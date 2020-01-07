using namespace std;
using namespace boost;

// list of scheduled jobs
// list of running jobs

uint64_t current_vertex;
set<DirectedGraph::vertex_descriptor> active_nodes;
set<DirectedGraph::vertex_descriptor> running_nodes;
priority_queue<DirectedGraph::vertex_descriptor, vector<DirectedGraph::vertex_descriptor>, Compare> active_nodes_queue;
auto lb_total_t = omp_get_wtime();

struct children_set_bfs_visitor : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        G[current_vertex].children.insert(target(e, g));
    }
};

void init_children_set(DirectedGraph::vertex_descriptor vertex) {
    current_vertex = vertex;
    children_set_bfs_visitor vis;
    breadth_first_search(G, vertex, visitor(vis));
}

void lb_init_graph() {
#ifdef DEBUG_PRINT
    cout << "Initializing vertex properties..." << endl;
#endif

    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    for (; vertexIt != vertexEnd; ++vertexIt) {
        G[*vertexIt].level = 0;
        G[*vertexIt].is_predicate_node = false;
        G[*vertexIt].is_active = false;
        G[*vertexIt].proc_time = 0;
    }
}

/*
 * These functions do not work yet but essentially they are
 * replicating LB construction of interval lists
 */
void set_interval_lists(DirectedGraph::vertex_descriptor vertex) {
    DirectedGraph::in_edge_iterator inedgeIt, inedgeEnd;
    tie(inedgeIt, inedgeEnd) = in_edges(vertex, G);
    for (; inedgeIt != inedgeEnd; ++inedgeIt) {
        G[source(*inedgeIt, G)].int_interval += G[vertex].int_interval;
        DirectedGraph::vertex_descriptor vert = source(*inedgeIt, G);
        set_interval_lists(vert);
    }
}

void leaf_int_set(size_t max_top_id) {
    auto count = 0;
    vector<DirectedGraph::vertex_iterator> leafs;
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd ; ++vertexIt) {
        if (G[*vertexIt].level == 0) {
#ifdef DEBUG_PRINT
            cout << "Found node with level 0. Skipping it..." << endl;
#endif
            continue;
        }
        G[*vertexIt].reach_id = max_top_id - G[*vertexIt].top_id;
        G[*vertexIt].int_interval.insert(G[*vertexIt].reach_id);
        DirectedGraph::out_edge_iterator outedgeIt, outedgeEnd;
        tie(outedgeIt, outedgeEnd) = out_edges(*vertexIt, G);
        auto num_outedges = 0;
        for (; outedgeIt != outedgeEnd; ++outedgeIt) {
            num_outedges++;
        }
        if (num_outedges == 0) {
            leafs.push_back(vertexIt);
            count++;
        }
    }
#ifdef DEBUG_PRINT
    cout << "Number of leafs: " << leafs.size() << endl;
#endif
    for (auto vertex : leafs) {
        set_interval_lists(*vertex);
    }
#ifdef DEBUG_PRINT
    cout << "Done setting interval lists." << endl;
#endif
}

// TODO: this function is broken. Fix it.
void construct_interval_list() {
    vector<graph_traits<DirectedGraph>::vertex_descriptor> cont;
    auto start_time = omp_get_wtime();
    try {
        topological_sort(G, std::back_inserter(cont));
    } catch (not_a_dag) {
        cout << "Graph is not a DAG." << endl;
        return ;
    }
    auto end_time = omp_get_wtime();
#ifdef DEBUG_PRINT
    cout << "Topological sort took: " << (end_time - start_time) << " seconds." << endl;
#endif

    cout << "A topological ordering:" << endl;
    size_t j = 0;
    for (vector<graph_traits<DirectedGraph>::vertex_descriptor>::reverse_iterator ii = cont.rbegin() + 1; ii != cont.rend(); ++ii) {
        G[*ii].top_id = j++;
    }
    cout << "Setting up children set." << endl;
    start_time = omp_get_wtime();
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if (G[*vertexIt].level == 0) {
            continue;
        }
        init_children_set(*vertexIt);
    }
    cout << "Done with children set." << endl;

    tie(vertexIt, vertexEnd) = vertices(G);
    return ;

//    leaf_int_set(j - 1);
    end_time = omp_get_wtime();
#ifdef DEBUG_PRINT
    cout << "Constructing interval lists took: " << (end_time - start_time) << " seconds." << endl;
#endif
}
/* Functions below are new and supposed to work */

void activate_lb_children(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
    tie(edgeIt, edgeEnd) = out_edges(node, G);
    for (; edgeIt != edgeEnd; ++edgeIt) {
        if (G[target(*edgeIt, G)].is_active) {
            DirectedGraph::out_edge_iterator outedgeIt, outedgeEnd;
            tie(outedgeIt, outedgeEnd) = out_edges(target(*edgeIt, G), G);
            for (; outedgeIt != outedgeEnd; ++outedgeIt) {
                active_nodes.insert(target(*outedgeIt, G));
            }
        }
    }
}

DirectedGraph::vertex_descriptor find_node() {
    for (auto anode : active_nodes) {
        int found = 0;
        // Search among running nodes
        for (auto rnode : running_nodes) {
            if (G[rnode].children.find(anode) != G[rnode].children.end()) {
                found = 1;
                break;
            }
        }
        // Search among active nodes
        if (!found) {
            for (auto rnode : active_nodes) {
                if (G[rnode].children.find(anode) != G[rnode].children.end()) {
                    found = 1;
                    break;
                }
            }
        }
        if (!found) {
            return anode;
        }
    }
}

void scheduler(double &counter) {
    auto start_t = omp_get_wtime();
    while(!active_nodes.empty() && (running_nodes.size() <= 8)) {
        set<DirectedGraph::vertex_descriptor>::iterator it = active_nodes.begin();
        DirectedGraph::vertex_descriptor node = find_node();
        running_nodes.insert(node);
        active_nodes.erase(node);
        counter++;
    }
    auto end_t = omp_get_wtime();
    lb_total_t += (end_t - start_t);
}

// TODO:
// create a class that will represent a working thread such that
// it will hold current job running with its time stamp. Then scheduler
// can go through the list of current running jobs to find, which next
// job can be scheduled and then replace job with earliest end time and
// update running job/time stamp on that thread accordingly. While updating
// also need to update set of active jobs and be careful not to run same job
// twice as priority queue is a vector but not set.

void run_lb_scheduler() {
    double counter = 0;
#ifdef DEBUG_PRINT
    cout << "Starting LB scheduler." << endl;
#endif
    lb_total_t = 0;
    scheduler(counter);
    while(!active_nodes.empty() || !running_nodes.empty()) {
        if (!running_nodes.empty()) {
            set<DirectedGraph::vertex_descriptor>::iterator it = running_nodes.begin();
            activate_lb_children(*it);
            running_nodes.erase(*it);
            scheduler(counter);
        } else {
            cout << "BUG: There should be nodes running." << endl;
            exit(1);
        }
    }
#ifdef DEBUG_PRINT
    cout << "LB Scheduler took " << lb_total_t << " seconds." << endl;
    cout << "Approximate time in units to run " << counter << " jobs is: " << ceil(counter / num_proc) << endl;
#endif
}

void lb_populate_active_list_from_trace() {
    ifstream start_node;
    start_node.open("startingNodes.dat");
    if (!start_node) {
        cout << "Error during opening startingNodes.dat file." << endl;
        exit(ENOENT);
    }

    size_t num_lines = count_lines_in_file("startingNodes.dat");
    string node_orig;
    for (auto i = 0; i < num_lines; ++i) {
        start_node >> node_orig;
        DirectedGraph::vertex_descriptor vertex = vertices_mapping[node_orig];
        DirectedGraph::out_edge_iterator eit, eend;
        tie(eit, eend) = out_edges(vertex, G);
        for (; eit != eend; ++eit) {
            if (!G[target(*eit, G)].is_active) {
                cout << "One of the children of start node is not active." << endl;
                exit(1);
            } else {
                active_nodes_queue.push(target(*eit, G));
            }
        }
    }
#ifdef DEBUG_PRINT
    cout << "All children of start nodes are active." << endl;
#endif
}

void lb_populate_active_list() {
    auto count = 0;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        auto num_in_edges = 0;
        DirectedGraph::in_edge_iterator eit, eend;
        for (; eit != eend; ++eit) {
            num_in_edges++;
        }
        if (num_in_edges > 0) {
            continue;
        }

        DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
        tie(edgeIt, edgeEnd) = out_edges(*vit, G);
        for (; edgeIt != edgeEnd; ++edgeIt) {
            if (G[target(*edgeIt, G)].is_active) {
                DirectedGraph::out_edge_iterator outedgeIt, outedgeEnd;
                tie(outedgeIt, outedgeEnd) = out_edges(target(*edgeIt, G), G);
                for (; outedgeIt != outedgeEnd; ++outedgeIt) {
                    active_nodes.insert(target(*outedgeIt, G));
                    count++;
                }
            }
        }
    }
#ifdef DEBUG_PRINT
    cout << "Added " << count << " jobs to initial active list." << endl;
#endif
}

void lb_scheduler(bool trace_enabled) {
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    lb_init_graph();
    auto start_t = omp_get_wtime();
    for (; vit != vend; ++vit) {
        init_children_set(*vit);
    }
    auto end_t = omp_get_wtime();
#ifdef DEBUG_PRINT
    cout << "Constructing children set took " << (end_t - start_t) << " seconds." << endl;
#endif
    if (trace_enabled) {
#ifdef DEBUG_PRINT
        check_trace();
        size_t num_active_jobs = 0;
        size_t num_elements = 0;
        tie(vit, vend) = vertices(G);
        for (; vit != vend; ++vit) {
            num_elements++;
            if (G[*vit].is_active) {
                num_active_jobs++;
            }
        }
        cout << "Number of active jobs (check): " << num_active_jobs << endl;
        cout << "Number of elements (check): " << num_elements << endl;
#endif
        assign_jobs_from_trace();
    } else {
        assign_jobs();
    }
    if (trace_enabled) {
        lb_populate_active_list_from_trace();
    } else {
        lb_populate_active_list();
    }

    run_lb_scheduler();
}
