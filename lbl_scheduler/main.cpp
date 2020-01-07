#include "defs.hpp"
#include "sbu_scheduling.hpp"
#include "lb_scheduling.hpp"
#include "stashed_code.hpp"

int main(int argc, char **argv) {
    ifstream edges_file;
    edges_file.open("edges.dat");

    if(!edges_file) {
        cout << "Error during opening the file." << endl;
        return ENOENT;
    }

    size_t num_edges_in_file = count_lines_in_file("edges.dat");
    string orig_src, orig_dst;
    uint64_t src, dst;
    vector<DirectedGraph::vertex_descriptor> pred_nodes;
    int c = 0;
    uint64_t num_init_tasks = 100;
    bool trace_enabled = false;
    
    while((c = getopt(argc, argv, "l:g:p:t")) != -1) {
        switch(c) {
            case 'g':
                num_init_tasks = atoi(optarg);
                break;
            case 'p':
                num_proc = atoi(optarg);
                break;
            case 't':
                trace_enabled = true;
                break;
            case 'l':
                lookahead_levels = atoi(optarg);
                break;
            default:
                cout << "Error during parsing flags." << endl;
                return EINVAL;
        }
    }
#ifdef DEBUG_PRINT
    cout << "Starting parsing the file." << endl;
#endif
    auto start_time = omp_get_wtime();
    for (auto i = 0; i < num_edges_in_file; ++i) {
        edges_file >> orig_src >> orig_dst >> src >> dst;
        vertices_mapping[orig_src] = src;
        vertices_mapping[orig_dst] = dst;
        reverse_vertices_mapping[to_string(src)] = orig_src;
        reverse_vertices_mapping[to_string(dst)] = orig_dst;
        add_edge(src, dst, G);
    }

    // Setting initial values for most of the vector properties
    init_graph();

    auto end_time = omp_get_wtime();
    edges_file.close();
#ifdef DEBUG_PRINT
    cout << "Graph creation finished in: " << (end_time - start_time) << " seconds." << endl;
    cout << "Number of edges: " << num_edges(G) << endl;
    cout << "Number of vertices: " << num_vertices(G) << endl;
    cout << "Number of leafs in the graph: " << count_leafs() << endl;
    cout << "Number of processors of simulated system: " << num_proc << endl;
#endif

    if (!is_dag(G)) {
        cout << "Graph is not a DAG." << endl;
        return EINVAL;
    }
#ifdef DEBUG_PRINT
    cout << "Graph is a DAG." << endl;
#endif

    ifstream prednodes_file;
    prednodes_file.open("prednodes.dat");

    if(!prednodes_file) {
        cout << "Error during opening the file with list of predicate nodes." << endl;
        return ENOENT;
    }

    num_pred_nodes = count_lines_in_file("prednodes.dat");
#ifdef DEBUG_PRINT
    cout << "Number of predicate nodes: " << num_pred_nodes << endl;
#endif
    for (auto i = 0; i < num_pred_nodes; ++i) {
        prednodes_file >> orig_src;
        G[vertices_mapping[orig_src]].is_predicate_node = true;
        pred_nodes.push_back(vertices_mapping[orig_src]);
    }
    prednodes_file.close();

#ifdef ACTIVATE
    vector<DirectedGraph::vertex_descriptor> pred_nodes_source;
    for (auto node : pred_nodes) {
        auto num_in_edges = 0;
        DirectedGraph::in_edge_iterator eit, eend;
        tie(eit, eend) = in_edges(node, G);
        for (; eit != eend; ++eit) {
            num_in_edges++;
        }
        if (num_in_edges == 0) {
            pred_nodes_source.push_back(node);
            G[node].is_source_node = true;
        }
        if (num_in_edges > 1) {
            cout << "More than one incoming edge for predicate node " << node << endl;
            // TODO: currently just exit the program.
            exit(1);
        }
    }

    vector<DirectedGraph::vertex_descriptor> unit_nodes_source;
    vector<DirectedGraph::vertex_descriptor> nodes_source;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        auto num_in_edges = 0;
        DirectedGraph::in_edge_iterator eit, eend;
        tie(eit, eend) = in_edges(*vit, G);
        for (; eit != eend; ++eit) {
            num_in_edges++;
        }
        if (num_in_edges == 0) {
            if (!G[*vit].is_predicate_node) {
                unit_nodes_source.push_back(*vit);
                cout << "WARNING: Source unit node." << endl;
            }
            nodes_source.push_back(*vit);
            G[*vit].is_source_node = true;
        }
    }
#ifdef DEBUG_PRINT
    cout << "All predicate nodes have at most one incoming edge." << endl;
#endif
    check_leafs(reverse_vertices_mapping);
    #ifdef DEBUG_PRINT
    cout << "Number of source nodes: " << nodes_source.size() << endl;
    cout << "Number of source predicate nodes: " << pred_nodes_source.size() << endl;
    cout << "Number of source unit nodes: " << unit_nodes_source.size() << endl;
    #endif
    if (random_activation(num_init_tasks, pred_nodes_source) < 0) {
    #ifdef DEBUG_PRINT
        cout << "Error during random activation." << endl;
    #endif
        return EINVAL;
    }
    if (save_active_graph(reverse_vertices_mapping) < 0) {
    #ifdef DEBUG_PRINT
        cout << "Error during saving activated graph into the file." << endl;
    #endif
        return ENOENT;
    }
    #ifdef DEBUG_PRINT
    auto num_active_nodes = 0;
    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        if (G[*vit].is_active) {
            num_active_nodes++;
        }
    }
    cout << "Number of activated nodes: " << num_active_nodes << endl;
    cout << "Random activation complete." << endl;
    #endif
    return 0;
#endif

    /* Code above concludes initial graph creation */

#ifdef DRAW
    // TODO: this function will not compile. Fix it.
    auto levels_to_draw = num_levels;
    if (draw_graph(dummy_vertex, levels, levels_to_draw, "graph.gv") < 0) {
        cout << "Error during producing file for graph plotting." << endl;
        return -1;
    }
    #ifdef DEBUG_PRINT
    cout << "Drawing graph is complete. Please use dot to produce SVG file." << endl;
    #endif
#endif

    /* Place for a temporary code */
    /* Place for a temporary code */

#ifdef SBU_SCHEDULER
    sbu_scheduler(trace_enabled);
#elif LB_SCHEDULER
    lb_scheduler(trace_enabled);
#else
    cout << "Scheduler should be chosen." << endl;
    return -1;
#endif

    return 0;
}
