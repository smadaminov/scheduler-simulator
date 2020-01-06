#include <iostream>
#include <fstream>
#include <string>
#include <errno.h>
#include <omp.h>
#include <stdint.h>
#include <thread>
#include <mutex>
#include <chrono>
#include <queue>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/subgraph.hpp>

#include <boost/icl/interval_set.hpp>
#include <boost/icl/closed_interval.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/icl/left_open_interval.hpp>
#include <boost/icl/open_interval.hpp>

using namespace std;
using namespace boost;

struct vertex_properties {
    uint64_t level;
    uint64_t level_mod;
    bool is_predicate_node;
    bool is_active;
    bool is_scheduled;
    icl::interval_set<uint64_t> int_interval;
    uint64_t top_id;
    uint64_t reach_id;
    set<uint64_t> children;
    bool set_for_lb;
    double proc_time;
    bool is_finished;
    bool is_source_node;
    uint64_t in_degree;
    uint64_t out_degree;
    uint64_t seq_number;
    bool prunned;
    bool in_busy_watch;
    bool is_executed;
};

typedef adjacency_list<vecS, vecS, bidirectionalS, vertex_properties, no_property> DirectedGraph;

DirectedGraph G, G_tmp;
int activation_probability = 50;
int num_proc = 24;
int lookahead_levels = 0;
bool print_executed_nodes = false;
bool parse_ancestor_list = false;
mutex jobs_mutex;

map<string, uint64_t> vertices_mapping;
map<string, string> reverse_vertices_mapping;
vector<vector<DirectedGraph::vertex_iterator>> levels;

auto cmp = [](DirectedGraph::vertex_descriptor a, DirectedGraph::vertex_descriptor b) {
    return G[a].seq_number < G[b].seq_number;
};

// TODO: having scheduled_nodes set with the seq_number based comparator may influence
// the sbi scheduler. check that
set<DirectedGraph::vertex_descriptor, decltype(cmp)> scheduled_nodes(cmp);
set<DirectedGraph::vertex_descriptor, decltype(cmp)> running_nodes(cmp);
set<DirectedGraph::vertex_descriptor, decltype(cmp)> run_queue(cmp);
vector<set<DirectedGraph::vertex_descriptor>> ancestral_list;
vector<set<DirectedGraph::vertex_descriptor>> active_nodes_on_level;

size_t num_pred_nodes;
uint64_t num_executed_jobs;
uint64_t sink_node;

uint64_t current_level;
uint64_t num_run_nodes_on_level;

auto sched_start = chrono::high_resolution_clock::now();
auto sched_end = chrono::high_resolution_clock::now();
auto scheduling_overhead = chrono::duration_cast<chrono::microseconds>(sched_start - sched_start);

struct activating_bfs_visitor : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        if (G[source(e, G)].is_active) {
            auto r = ((double)rand()) / RAND_MAX;
            if (r < activation_probability) {
                G[target(e, G)].is_active = true;
            }
        }
    }
};

class Compare {
    public:
        bool operator() (DirectedGraph::vertex_descriptor& x, DirectedGraph::vertex_descriptor& y) {
            int comp = memcmp(&G[x].seq_number, &G[y].seq_number, sizeof(uint64_t));
            return comp > 0;
        }
};
            
size_t count_lines_in_file(const string &filename) {
    size_t number_of_lines = 0;
    ifstream edges_file;
    edges_file.open(filename);

    if(!edges_file) {
        cout << "Error during opening edges.dat file." << endl;
        return ENOENT;
    }

    string tmp_string;
    if (edges_file.is_open()) {
        istream& in = edges_file;
        while(std::getline(in, tmp_string)) {
            number_of_lines++;
        }
    } else {
        cout << "Attempt to count lines in the file failed. File is not open." << endl;
        return -ENOENT;
    }

    return number_of_lines;
}

void init_graph() {
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    for (; vertexIt != vertexEnd; ++vertexIt) {
        G[*vertexIt].level = -1;
        G[*vertexIt].level_mod = -1;
        G[*vertexIt].is_predicate_node = false;
        G[*vertexIt].is_active = false;
        G[*vertexIt].is_scheduled = false;
        G[*vertexIt].set_for_lb = false;
        G[*vertexIt].proc_time = 0;
        G[*vertexIt].is_finished = false;
        G[*vertexIt].is_source_node = false;
        G[*vertexIt].in_degree = -1;
        G[*vertexIt].out_degree = -1;
        G[*vertexIt].seq_number = -1;
        G[*vertexIt].in_busy_watch = false;
        G[*vertexIt].is_executed = false;
    }
}

void set_node_priorities() {
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    ifstream sequence_file;
    sequence_file.open("sequence.dat");

    uint64_t orig_node_id, node_id, priority;
    size_t num_lines_in_file = count_lines_in_file("sequence.dat");
    for (auto i = 0; i < num_lines_in_file; ++i) {
        sequence_file >> orig_node_id >> node_id >> priority;
        G[node_id].seq_number = priority;
    }

    sequence_file.close();
}

bool is_dag(DirectedGraph G) {
    vector<graph_traits<DirectedGraph>::vertex_descriptor> cont;
    try {
        topological_sort(G, std::back_inserter(cont));
    } catch (not_a_dag) {
        return false;
    }

    return true;
}

uint64_t find_depth() {
    uint64_t max_level = 0;
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if (G[*vertexIt].level > max_level) {
            max_level = G[*vertexIt].level;
        }
    }

    return max_level;
}

int draw_graph(uint64_t dummy_vertex, const vector<vector<DirectedGraph::vertex_iterator> > &levels, uint64_t levels_to_draw, const string &output_filename) {
    cout << "Drawing graph with " << levels_to_draw << " levels." << endl;
    ofstream output_file;
    output_file.open(output_filename);

    if (!output_file) {
        cout << "Error during opening the file for outputing graph in the DOT/Graphviz format." << endl;
        return -ENOENT;
    }

    string indent = "  ";
    output_file << "digraph G {" << endl
                << indent << "rankdir = TB;" << endl
                << indent << "graph [nodesep=\"1\", ranksep=\"5\"];" << endl
                << indent << "splines=\"false\";" << endl
                << indent << "node[shape = square];" << endl
                << indent << "subgraph {" << endl;

    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    DirectedGraph::out_edge_iterator outedgeIt, outedgeEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    auto count_nodes = 0;
    auto count_edges = 0;
    auto skipped = 0;
    cout << "There are " << num_vertices(G) << " nodes." << endl;
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if (G[*vertexIt].prunned) {
            skipped++;
            continue;
        }
        count_nodes++;
#ifdef DEBUG_PRINT
//        cout << "Entering edges from: " << *vertexIt << endl;
#endif
        tie(outedgeIt, outedgeEnd) = out_edges(*vertexIt, G);
        for (; outedgeIt != outedgeEnd; ++outedgeIt) {
            if (G[target(*outedgeIt, G)].level <= levels_to_draw) {
                output_file << indent << indent << *vertexIt << " -> " << target(*outedgeIt, G) << endl;
                count_edges++;
            }
        }
    }

#ifdef DEBUG_PRINT
    cout << "Edges are entered." << endl;
    cout << "There are " << count_nodes << " nodes." << endl;
    cout << "There are " << count_edges << " edges." << endl;
    cout << "Skipped " << skipped << " nodes." << endl;
#endif

    for (auto i = 1; i < levels_to_draw; ++i) {
        if (levels[i].size() != 0) {
            output_file << indent << indent << "{ rank = same; ";
            for (auto vi : levels[i]) {
                if (G[*vi].prunned) {
                    cout << "ERROR: pruned node got level." << endl;
                    exit(1);
                }
                output_file << *vi << "; ";
            }
            output_file << "}" << endl;
        }
    }

#ifdef DEBUG_PRINT
    cout << "Levels (ranks) are entered." << endl;
#endif

    output_file << indent << "}" << endl << "}" << endl;

    output_file.close();

    return 0;
}

size_t count_leafs() {
    auto count = 0;
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd ; ++vertexIt) {
        if (out_degree(*vertexIt, G) == 0) {
            count++;
        }
    }

    return count;
}

void check_leafs(map<string, string> &reverse_vertices_mapping) {
#ifdef DEBUG_PRINT
    cout << "Checking leafs of the Graph." << endl;
#endif
    auto pred_count = 0;
    auto unit_count = 0;
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd ; ++vertexIt) {
        DirectedGraph::out_edge_iterator outedgeIt, outedgeEnd;
        tie(outedgeIt, outedgeEnd) = out_edges(*vertexIt, G);
        auto num_outedges = 0;
        for (; outedgeIt != outedgeEnd; ++outedgeIt) {
            num_outedges++;
        }
        if (num_outedges == 0) {
            if (G[*vertexIt].is_predicate_node) {
                pred_count++;
            } else {
                unit_count++;
            }
        }
    }
#ifdef DEBUG_PRINT
    if ((pred_count == 0) && (unit_count != 0)) {
        cout << "All leafs are unit nodes." << endl;
    } else if ((pred_count != 0) && (unit_count == 0)) {
        cout << "All leafs are predicate nodes." << endl;
    } else if ((pred_count != 0) && (unit_count != 0)) {
        cout << "Leafs are mixture of predicate and unit nodes." << endl;
    } else {
        cout << "Something went wrong during checking leafs." << endl;
        exit(1);
    }
    cout << "Finished checking leafs." << endl;
#endif
}

void generate_jobs(int num_init_tasks, vector<DirectedGraph::vertex_iterator> nodes, map<string, string> &reverse_vertices_mapping) {
    if (num_init_tasks <= 0) {
        cout << "Number of initial tasks should be greater than zero." << endl;
        return ;
    } else if (num_init_tasks > nodes.size()) {
        cout << "Number of initial tasksk should be less than number of source nodes." << endl;
        return ;
    }

    random_shuffle(nodes.begin(), nodes.end());
    for (auto i = 0; i < num_init_tasks; ++i) {
        // TODO: generating/activating should be algo independent
//        active_nodes_on_level[1].insert(*nodes[i]);
        G[*nodes[i]].is_active = true;
        activating_bfs_visitor vis;
        breadth_first_search(G, *nodes[i], visitor(vis));
    }

    ofstream output_file;
    output_file.open("active_nodes.dat");

    if (!output_file) {
        cout << "Error during opening the file for outputing active nodes." << endl;
        return ;
    }

    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    auto num_active_nodes = 0;
    auto num_inactive_nodes = 0;
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if (G[*vertexIt].is_active) {
            num_active_nodes++;
            output_file << *vertexIt << " " << reverse_vertices_mapping[to_string(*vertexIt)] << endl;
        } else {
            num_inactive_nodes++;
        }
    }

#ifdef DEBUG_PRINT
    cout << "Number of active nodes: " << num_active_nodes << endl;
    cout << "Number of inactive nodes: " << num_inactive_nodes << endl;
#endif

    output_file.close();
}

void add_start_nodes() {
    ifstream input_file;
    input_file.open("startingNodes.dat");

    uint64_t node;
    string node_orig;

    auto num_active_nodes = 0;
    size_t num_lines_in_file = count_lines_in_file("startingNodes.dat");

    for (auto i = 0; i < num_lines_in_file; ++i) {
        input_file >> node_orig;
        if (G[vertices_mapping[node_orig]].level != 0) {
            cout << "ERROR: starting node is not source node." << endl;
            exit(1);
        }
        G[vertices_mapping[node_orig]].is_active = true;
        num_active_nodes++;
    }

#ifdef DEBUG_PRINT
    cout << "Number of starting nodes (from file): " << num_active_nodes << endl;
#endif

    input_file.close();
}

void assign_jobs_from_trace() {
    ifstream input_file;
    input_file.open("activatedJob.dat");

    uint64_t node;
    string node_orig;

    auto num_active_nodes = 0;
    size_t num_lines_in_file = count_lines_in_file("activatedJob.dat");

    for (auto i = 0; i < num_lines_in_file; ++i) {
        input_file >> node_orig;
        map<string, uint64_t>::iterator it = vertices_mapping.find(node_orig);
        if (it == vertices_mapping.end()) {
            cout << "ERROR: Key not found (during setting activated jobs): " << node_orig << endl;
            exit(1);
        }
        G[vertices_mapping[node_orig]].is_active = true;
        if (in_degree(vertices_mapping[node_orig], G) == 0) {
            cout << "ERROR: Activating source node." << endl;
            exit(1);
        }
        num_active_nodes++;
    }

#ifdef DEBUG_PRINT
    cout << "Number of active nodes (from file): " << num_active_nodes << endl;
#endif

    input_file.close();

    double proc_time;
    input_file.open("processingTime.dat");
    num_lines_in_file = count_lines_in_file("processingTime.dat");

    for (auto i = 0; i < num_lines_in_file; ++i) {
        input_file >> node_orig >> proc_time;
        map<string, uint64_t>::iterator it = vertices_mapping.find(node_orig);
        if (it == vertices_mapping.end()) {
            cout << "ERROR: Key not found (during setting processing times): " << node_orig << endl;
            exit(1);
        }
        G[vertices_mapping[node_orig]].proc_time = proc_time;
    }

    input_file.close();

#ifdef DEBUG_PRINT
    cout << "Jobs have been assigned." << endl;
#endif
}

void assign_jobs() {
    ifstream input_file;
    input_file.open("active_nodes.dat");

    uint64_t node;
    string node_orig;
    bool set_for_lb;

    auto num_active_nodes = 0;
    size_t num_lines_in_file = count_lines_in_file("active_nodes.dat");

    for (auto i = 0; i < num_lines_in_file; ++i) {
        input_file >> node >> node_orig >> set_for_lb;
        G[node].is_active = true;
        G[node].set_for_lb = set_for_lb;
        num_active_nodes++;
    }

#ifdef DEBUG_PRINT
    cout << "Number of active nodes (from file): " << num_active_nodes << endl;
#endif

    input_file.close();
}

int save_active_graph(map<string, string> &reverse_vertices_mapping) {
    ofstream output_file;
    output_file.open("active_nodes.dat");

    if (!output_file) {
        cout << "Error during opening the file for outputing active nodes." << endl;
        return -ENOENT;
    }

    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    auto num_active_nodes = 0;
    auto num_inactive_nodes = 0;
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if (G[*vertexIt].is_active) {
            num_active_nodes++;
            output_file << *vertexIt << " " << reverse_vertices_mapping[to_string(*vertexIt)] << " " << G[*vertexIt].set_for_lb << endl;
        } else {
            num_inactive_nodes++;
        }
    }

#ifdef DEBUG_PRINT
    cout << "Number of active nodes: " << num_active_nodes << endl;
    cout << "Number of inactive nodes: " << num_inactive_nodes << endl;
#endif

    output_file.close();

    return 0;
}

void activate_children(uint64_t node) {
    DirectedGraph::out_edge_iterator outedgeIt, outedgeEnd;
    tie(outedgeIt, outedgeEnd) = out_edges(node, G);
    for (; outedgeIt != outedgeEnd; ++outedgeIt) {
        auto r = rand() % 100;
        if (G[node].is_predicate_node) {
            if (G[node].is_active) {
                r = 100;
            } else {
                r = 0;
            }
        }
        if (!G[node].is_predicate_node) {
            if (!G[node].is_active) {
                r = 0;
            }
        }
        if (r > activation_probability) {
            G[target(*outedgeIt, G)].is_active = true;
            activate_children(target(*outedgeIt, G));
        }
    }
}

int random_activation(uint64_t num_init_tasks, vector<DirectedGraph::vertex_descriptor> nodes) {
    if (num_init_tasks <= 0) {
#ifdef DEBUG_PRINT
        cout << "Number of initial tasks should be greater than zero." << endl;
#endif
        return -1;
    } else if (num_init_tasks > nodes.size()) {
#ifdef DEBUG_PRINT
        cout << "Number of initial tasks should be less than number of predicate nodes." << endl;
#endif
        return -2;
    }

    srand(1);
    random_shuffle(nodes.begin(), nodes.end());
    for (auto i = 0; i < num_init_tasks; ++i) {
        G[nodes[i]].is_active = true;
        G[nodes[i]].set_for_lb = true;
        activate_children(nodes[i]);
    }

    return 0;
}

void rebuild_graph() {
    cout << "Rebuilding graph..." << endl;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    size_t num_pruned = 0;
    for ( ; vit != vend; ++vit) {
        if (G[*vit].is_source_node) {
            continue;
        }
        if (G[*vit].is_predicate_node) {
            DirectedGraph::in_edge_iterator in_edge_it, in_edge_end;
            tie(in_edge_it, in_edge_end) = in_edges(*vit, G);
            DirectedGraph::out_edge_iterator eit, eend;
            tie(eit, eend) = out_edges(*vit, G);
            for ( ; eit != eend; ++eit) {
                add_edge(source(*in_edge_it, G), target(*eit, G), G);
            }
            remove_edge(source(*in_edge_it, G), *vit, G);
            clear_out_edges(*vit, G);
            if (G[*vit].is_active) {
                cout << "ERROR: Active non-source predicate node." << endl;
                exit(1);
            }
            num_pruned++;
            G[*vit].prunned = true;
        }
    }

    cout << "Graph has been rebuilt." << endl;
    cout << "Number of vertices pruned: " << num_pruned << endl;
    cout << "Number of edges: " << num_edges(G) << endl;
    cout << "Number of vertices: " << num_vertices(G) - num_pruned << endl;
    copy_graph(G, G_tmp);
}

void add_sink_node() {
    add_vertex(G);
    sink_node = num_vertices(G) - 1;
    G[sink_node].proc_time = 0;
    G[sink_node].level = 0;
    G[sink_node].is_scheduled = false;
}

uint64_t find_min_running_level() {
    uint64_t min_run_level = numeric_limits<uint64_t>::max();

    for (auto node : running_nodes) {
        if (G[node].level < min_run_level)
            min_run_level = G[node].level;
    }

    return min_run_level;
}

/*
void elapsed_time(auto start, auto end) {
    cout << "Elapsed time in microseconds : " 
        << chrono::duration_cast<chrono::microseconds>(end - start).count()
        << " µs" << endl;

    cout << "Elapsed time in milliseconds : " 
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << endl;

    cout << "Elapsed time in seconds : " 
        << chrono::duration_cast<chrono::seconds>(end - start).count()
        << " sec";
}
*/
