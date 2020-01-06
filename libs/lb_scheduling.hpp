using namespace std;
using namespace boost;

// list of scheduled jobs
// list of running jobs

uint64_t current_vertex;
uint64_t num_start_nodes;
set<DirectedGraph::vertex_descriptor> active_nodes;
set<DirectedGraph::vertex_descriptor> start_nodes;

priority_queue<DirectedGraph::vertex_descriptor, vector<DirectedGraph::vertex_descriptor>, Compare> active_nodes_queue;
auto lb_total_t = omp_get_wtime();

DirectedGraph::vertex_descriptor lookup_node;
bool lookup_node_found;
bool ancestor_list_converged = false;

set<DirectedGraph::vertex_descriptor, decltype(cmp)> activation_queue(cmp);
set<DirectedGraph::vertex_descriptor> busy_set;
//vector<set<DirectedGraph::vertex_descriptor, decltype(cmp)>> busy_watch_set;
vector<set<DirectedGraph::vertex_descriptor>> busy_watch_set;

double current_timestamp;

ofstream executed_nodes_file;

struct ancestral_bfs : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        if (target(e, G) == lookup_node) {
            lookup_node_found = true;
        }
    }
};

struct merge_ancestral_info_bfs : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        auto node = target(e, G);
        DirectedGraph::in_edge_iterator inedgeIt, inedgeEnd;
        tie(inedgeIt, inedgeEnd) = in_edges(node, G);

        auto list_size = ancestral_list[node].size();

        for (; inedgeIt != inedgeEnd; ++inedgeIt) {
            // TODO: Re-write using std::merge to make the insertion more efficient
            ancestral_list[node].insert(ancestral_list[source(*inedgeIt, G)].begin(), ancestral_list[source(*inedgeIt, G)].end());
            ancestral_list[node].insert(source(*inedgeIt, G));
        }

        if (list_size != ancestral_list[node].size()) {
            ancestor_list_converged = false;
        }
    }
};

struct children_set_bfs_visitor : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        G[current_vertex].children.insert(target(e, g));
    }
};

void finish_job_h(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree, int cpu_num);

void init_children_set(DirectedGraph::vertex_descriptor vertex) {
    current_vertex = vertex;
    children_set_bfs_visitor vis;
    breadth_first_search(G, vertex, visitor(vis));
}

void lb_init_graph(vector<DirectedGraph::vertex_iterator> &source_nodes) {
#ifdef DEBUG_PRINT
    cout << "Initializing vertex properties and constructing set of source nodes..." << endl;
#endif
    uint64_t source_unit_nodes = 0;
    num_start_nodes = 0;
    num_executed_jobs = 0;

    if (print_executed_nodes)
        executed_nodes_file.open("executed_nodes.dat");

    DirectedGraph::in_edge_iterator inedgeIt, inedgeEnd;
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    current_level = 1;
    num_run_nodes_on_level = 0;

    for (; vertexIt != vertexEnd; ++vertexIt) {
        G[*vertexIt].level = -1;
        G[*vertexIt].proc_time = 0;
        G[*vertexIt].level_mod = -1;
        G[*vertexIt].is_active = false;
        G[*vertexIt].is_scheduled = false;
        G[*vertexIt].in_degree = in_degree(*vertexIt, G);
        G[*vertexIt].prunned = false;
        if (G[*vertexIt].in_degree == 0) {
            source_nodes.push_back(vertexIt);
            G[*vertexIt].is_source_node = true;
            if (!G[*vertexIt].is_predicate_node) {
                cout << "WARNING: Source unit node.\r";
//                cout << "WARNING: Source unit node: " << *vertexIt << endl;
                source_unit_nodes++;
            }
        }
    }

    if (source_unit_nodes > 0) {
        cout << endl;
    }

#ifdef DEBUG_PRINT
    cout << "Number of source nodes: " << source_nodes.size() << endl;
    cout << "Number of source unit nodes: " << source_unit_nodes << endl;
#endif
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

bool is_safe_to_run(DirectedGraph::vertex_descriptor node) {
    // Perform check among scheduled nodes
    for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
        if (*it == node)
            continue;
        if (ancestral_list[node].find(*it) != ancestral_list[node].end())
            return false;
    }
    // Perform check among running nodes
    for (auto it = running_nodes.begin(); it != running_nodes.end(); ++it) {
        if (*it == node) {
            cout << "ERROR: checked node is currently running." << endl;
            cout << "Node in question: " << node << endl;
            cout << "List of scheduled nodes:" << endl;
            for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
                cout << *it << endl;
            }
            cout << "List of running nodes:" << endl;
            for (auto it = running_nodes.begin(); it != running_nodes.end(); ++it) {
                cout << *it << endl;
            }
            exit(1);
        }
        if (ancestral_list[node].find(*it) != ancestral_list[node].end())
            return false;
    }

    return true;
}

bool is_safe_to_run_h(DirectedGraph::vertex_descriptor node) {
    // Perform check among scheduled nodes
    for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
        if (*it == node)
            continue;
        if (ancestral_list[node].find(*it) != ancestral_list[node].end())
            return false;
    }
    // Perform check among running nodes
    for (auto it = running_nodes.begin(); it != running_nodes.end(); ++it) {
        if (*it == node) {
            cout << "ERROR: checked node is currently running." << endl;
            cout << "Node in question: " << node << endl;
            cout << "List of scheduled nodes:" << endl;
            for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
                cout << *it << endl;
            }
            cout << "List of running nodes:" << endl;
            for (auto it = running_nodes.begin(); it != running_nodes.end(); ++it) {
                cout << *it << endl;
            }
            exit(1);
        }
        if (ancestral_list[node].find(*it) != ancestral_list[node].end())
            return false;
    }

    for (auto it = run_queue.begin(); it != run_queue.end(); ++it) {
        if (*it == node) {
            cout << "Node is already in the run queue." << endl;
            exit(1);
        }
        if (ancestral_list[node].find(*it) != ancestral_list[node].end()) {
            return false;
        }
    }

    return true;
}

DirectedGraph::vertex_descriptor get_scheduled_node() {
    sched_start = chrono::high_resolution_clock::now();
    DirectedGraph::vertex_descriptor node = -1;

//    for (auto i = 0; i < 1000; ++i) {
        for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
            if (is_safe_to_run(*it)) {
                node = *it;
                break;
            }
        }
//    }
    sched_end = chrono::high_resolution_clock::now();
    scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

    return node;
}

DirectedGraph::vertex_descriptor get_scheduled_node_h(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
again:
    sched_start = chrono::high_resolution_clock::now();
    DirectedGraph::vertex_descriptor node = -1;

    if (!run_queue.empty()) {
        node = *(run_queue.begin());

        run_queue.erase(node);

        sched_end = chrono::high_resolution_clock::now();
        scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

        return node;
    }
//    for (auto i = 0; i < 1000; ++i) {
        for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
            if (is_safe_to_run(*it)) {
                node = *it;
//                run_queue.insert(node);
//                scheduled_nodes.erase(node);
                break;
            }
        }

        for (auto it = run_queue.begin(); it != run_queue.end(); ++it) {
            scheduled_nodes.erase(*it);
        }

//    }
    sched_end = chrono::high_resolution_clock::now();
    scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

    /*
    auto sched_overhead = chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);
    double sched_overhead_c = sched_overhead.count();

    int counter_h = 0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] <= current_timestamp + sched_overhead_c / 1000.0 && run_on_proc[i] != sink_node) {
//            cout << "Finishing: " << run_on_proc[i] << endl;
            finish_job_h(end_time_on_proc, run_on_proc, CPUfree, allCPUfree, i);
            counter_h++;
        }
    }

    if (node != -1)
        return node;

    if (counter_h > 0) {
//        cout << "Checking again..." << endl;
        goto again;
    }
    */

    return node;
}

DirectedGraph::vertex_descriptor get_hybrid_scheduled_node() {
    sched_start = chrono::high_resolution_clock::now();
    DirectedGraph::vertex_descriptor node = -1;

    // First run nodes on the current level
    if (!run_queue.empty()) {
        node = *(run_queue.begin());
        run_queue.erase(node);

        /*
        // Here we perform check that everything
        // on the levels above is done
        int flag = 0;
        for (auto node : scheduled_nodes) {
            if (G[node].level < current_level) {
                cout << "ERROR: There are nodes on upper levels." << endl;
                cout << "Current level: " << current_level << endl;
                cout << "Node in question level: " << G[node].level << endl;
                cout << "Node: " << node << endl;
                flag = 1;
            }
        }

        if (flag) {
            cout << "List of scheduled nodes:" << endl;
            for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
                cout << *it << endl;
            }
            cout << "List of running nodes:" << endl;
            for (auto it = running_nodes.begin(); it != running_nodes.end(); ++it) {
                cout << *it << endl;
            }
//            exit(1);
        }
        */

        sched_end = chrono::high_resolution_clock::now();
        scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

        return node;
    }

///////////////////////////////////////////////
/*
    uint64_t tmp_level2 = current_level + 1;
    for (auto tmp_level = tmp_level2; tmp_level < tmp_level2 + lookahead_levels; tmp_level++) {
        if (tmp_level >= find_depth())
            return -1;
        for (auto l_node = active_nodes_on_level[tmp_level].begin(); l_node != active_nodes_on_level[tmp_level].end(); ++l_node) {
            if (scheduled_nodes.find(*l_node) == scheduled_nodes.end())
                continue;
            lookahead_bfs_visitor vis;
            next_node_found = false;
            next_node = *l_node;
            for (auto el : scheduled_nodes) {
                if (el == *l_node) {
                    continue;
                }
                breadth_first_search(G, el, visitor(vis));
            }
            if (next_node_found) {
                cout << "FOUND JOB ON LEVEL: " << tmp_level << endl;
                return *l_node;
            }
        }
    }

*/
//////////////////////////////////////////////

    // If there are no nodes on the current level then run LogicBlox scheduler
    for (auto it = scheduled_nodes.begin(); it != scheduled_nodes.end(); ++it) {
        if (is_safe_to_run(*it)) {
            node = *it;
            break;
        }
    }

    sched_end = chrono::high_resolution_clock::now();
    scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

    return node;
}

// TODO: we don't really need that function to be separate for the LB and SBU scheduler. Move it to the defs.hpp
// with some minor changes.
void lb_activate_children(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(node, G);
    for (; eit != eend; ++eit) {
        if (G[target(*eit, G)].is_active) {
            G[target(*eit, G)].is_scheduled = true;
            scheduled_nodes.insert(target(*eit, G));
            active_nodes_on_level[G[target(*eit, G)].level].insert(target(*eit, G));
        }
    }
}

void lb_activate_children_helper(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(node, G);
    // TODO: add a check against busy watch set
    for (; eit != eend; ++eit) {
        if (G[target(*eit, G)].is_active) {
            if (G[target(*eit, G)].is_predicate_node) {
                cout << "ERROR: Trying to activate predicate node." << endl;
                exit(1);
            }
            activation_queue.insert(target(*eit, G));
//            active_nodes_on_level[G[target(*eit, G)].level].insert(target(*eit, G));
        }
    }
}

void finish_earliest_job(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
    int cpu_num = -1;
    double min_time = numeric_limits<double>::max();
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            continue;
        if (end_time_on_proc[i] < min_time) {
            min_time = end_time_on_proc[i];
            cpu_num = i;
        }
    }

    if (cpu_num == -1) {
        cout << "ERROR: No job to finish." << endl;
        exit(1);
    }

    lb_activate_children(run_on_proc[cpu_num]);
    running_nodes.erase(run_on_proc[cpu_num]);

    // Once the active job finishes, the end time on all idle cpus
    // also should be updated
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            end_time_on_proc[i] = min_time;
    }

    current_timestamp = end_time_on_proc[cpu_num];

    if (print_executed_nodes)
        executed_nodes_file << run_on_proc[cpu_num] << " " << cpu_num << endl;

    run_on_proc[cpu_num] = sink_node;
    CPUfree = true;
    num_executed_jobs++;

    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] != sink_node)
            return ;
    }

    allCPUfree = true;
}

void finish_job_h(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree, int cpu_num) {
    if (cpu_num == -1) {
        cout << "ERROR: No job to finish." << endl;
        exit(1);
    }

    lb_activate_children(run_on_proc[cpu_num]);
    running_nodes.erase(run_on_proc[cpu_num]);

    // Once the active job finishes, the end time on all idle cpus
    // also should be updated
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            end_time_on_proc[i] = end_time_on_proc[cpu_num];
    }

    current_timestamp = end_time_on_proc[cpu_num];

    if (print_executed_nodes)
        executed_nodes_file << run_on_proc[cpu_num] << " " << cpu_num << endl;

    run_on_proc[cpu_num] = sink_node;
    CPUfree = true;
    num_executed_jobs++;

    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] != sink_node)
            return ;
    }

    allCPUfree = true;
}

void finish_hybrid_earliest_job(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
    cout << "Finishing earliest job..." << endl;
    int cpu_num = -1;
    double min_time = numeric_limits<double>::max();
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            continue;
        if (end_time_on_proc[i] < min_time) {
            min_time = end_time_on_proc[i];
            cpu_num = i;
        }
    }

    if (cpu_num == -1) {
        cout << "ERROR: No job to finish." << endl;
        exit(1);
    }

    lb_activate_children(run_on_proc[cpu_num]);
    running_nodes.erase(run_on_proc[cpu_num]);

    active_nodes_on_level[G[run_on_proc[cpu_num]].level].erase(run_on_proc[cpu_num]);

    if (G[run_on_proc[cpu_num]].level == current_level) {
        if (num_run_nodes_on_level > 1) {
            num_run_nodes_on_level--;
        } else {
            // TODO: optimize that min run level finding. It is better to
            // maintain a global variable to store that state
            uint64_t min_run_level = find_min_running_level();
            if (current_level == 27) {
                cout << "Moving from level 27 with size: " << active_nodes_on_level[current_level].size() << endl;
                cout << "Num run nodes: " << num_run_nodes_on_level << endl;
                cout << "Scheduled set size: " << scheduled_nodes.size() << endl;
            }
            if (min_run_level > current_level) {
                current_level++;
                if (!run_queue.empty()) {
                    cout << "ERROR: Run queue is not empty while moving to the next level." << endl;
//                    exit(1);
                }
                cout << "Going inside the level increment loop..." << endl;
                while(active_nodes_on_level[current_level].size() == 0) {
                    cout << "Empty level. Moving further..." << endl;
                    for (auto node : scheduled_nodes) {
                        if (G[node].level < current_level) {
                            cout << "ERROR: Level is supposed to be empty." << endl;
                            cout << "Current level: " << current_level << endl;
                            cout << "Node in question level: " << G[node].level << endl;
                        }
                    }
                    if (min_run_level <= current_level)
                        break;
                    current_level++;
                    if (current_level >= find_depth())
                        break;
                }
                cout << "Leaving the increment loop." << endl;
                cout << "Current level: " << current_level << endl;
                cout << "Max level: " << find_depth() << endl;
                if (current_level <= find_depth()) {
                    num_run_nodes_on_level = active_nodes_on_level[current_level].size();
                    for (auto node = active_nodes_on_level[current_level].begin(); node != active_nodes_on_level[current_level].end(); ++node)
                        cout << "NODE: " << *node << endl;
                    for (auto node : active_nodes_on_level[current_level]) {
                        cout << "Chimichanga: " << node << endl;
                        run_queue.insert(node);
                        cout << "Inserted: " << active_nodes_on_level[current_level].size() << endl;
                    }
                    cout << "Done with run queue insertion." << endl;
                }
            }
            cout << "Finished modifying current_level." << endl;
        }
    } else {
        uint64_t min_run_level = find_min_running_level();
        if (min_run_level > current_level) {
            current_level++;
            if (!run_queue.empty()) {
                cout << "ERROR: Run queue is not empty while moving to the next level." << endl;
//                exit(1);
            }
            cout << "Going inside the level increment loop... (2)" << endl;
            while(active_nodes_on_level[current_level].size() == 0) {
                cout << "Empty level. Moving further..." << endl;
                for (auto node : scheduled_nodes) {
                    if (G[node].level < current_level) {
                        cout << "ERROR: Level is supposed to be empty." << endl;
                        cout << "Current level: " << current_level << endl;
                        cout << "Node in question level: " << G[node].level << endl;
                    }
                }
                if (min_run_level <= current_level)
                    break;
                current_level++;
                if (current_level >= find_depth())
                    break;
            }
            cout << "Leaving the increment loop. (2)" << endl;
            if (current_level <= find_depth()) {
                num_run_nodes_on_level = active_nodes_on_level[current_level].size();
                for (auto node : active_nodes_on_level[current_level]) {
                    run_queue.insert(node);
                }
            }
        }
    }

    cout << "Done with leveling incremental." << endl;

    // Once the active job finishes, the end time on all idle cpus
    // also should be updated
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            end_time_on_proc[i] = min_time;
    }

    if (print_executed_nodes)
        executed_nodes_file << run_on_proc[cpu_num] << " " << cpu_num << endl;

    run_on_proc[cpu_num] = sink_node;
    CPUfree = true;
    num_executed_jobs++;

    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] != sink_node)
            return ;
    }
    allCPUfree = true;
}

void logic_blox_scheduler() {
    bool CPUfree = false;
    bool allCPUfree = false;
    vector<double> end_time_on_proc(num_proc, 0.0);
    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc, 0);

    if (scheduled_nodes.size() == 0) {
        cout << "No nodes scheduled." << endl;
        exit(1);
    }

    // Schedule first eight jobs on processors
    // TODO: fix a case when there is less than eight jobs initially available
    if (scheduled_nodes.size() < 8) {
        cout << "Less than eight initial nodes currently is not supported." << endl;
        exit(1);
    }

    // Scheduling first jobs on the CPUs
    for (auto i = 0; i < num_proc; ++i) {
//        auto node = get_scheduled_node();
        auto node = get_scheduled_node_h(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        if (node == -1) {
            for (auto j = i; j < num_proc; ++j) {
                end_time_on_proc[j] = 0;
                run_on_proc[j] = sink_node;
            }
            cout << "Beaker" << endl;
//            exit(1);
            break;
        }
        end_time_on_proc[i] = G[node].proc_time;
        run_on_proc[i] = node;
        scheduled_nodes.erase(node);
        running_nodes.insert(node);
    }

    /*
    current_time = numeric_limits<double>::max();
    for (auto i = 0; i < num_proc; ++i) {
        if (current_timestamp < end_time_on_proc[i])
            current_timestamp = end_time_on_proc[i];
    }
    */

    while (scheduled_nodes.size() != 0 || !allCPUfree) {
        if (!CPUfree) {
            finish_earliest_job(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        }
//        auto node = get_scheduled_node();
        auto node = get_scheduled_node_h(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        while (node == -1) {
            finish_earliest_job(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
//            node = get_scheduled_node();
            node = get_scheduled_node_h(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            if (allCPUfree && (node == -1))
                break;
        }
        if (allCPUfree && (node == -1))
            continue;

        int cpu_num = -1;
        double min_time = numeric_limits<double>::max();
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] != sink_node)
                continue;
            if (end_time_on_proc[i] < min_time) {
                min_time = end_time_on_proc[i];
                cpu_num = i;
            }
        }
        run_on_proc[cpu_num] = node;
        end_time_on_proc[cpu_num] += G[node].proc_time;
        scheduled_nodes.erase(node);
        running_nodes.insert(node);
        allCPUfree = false;
        CPUfree = false;
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] == sink_node) {
                CPUfree = true;
                break;
            }
        }
    }

    double max_time = -1.0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] > max_time) {
            max_time = end_time_on_proc[i];
        }
    }

    cout << "Total running time of schedule: " << max_time << endl;
    cout << "Number of executed jobs: " << num_executed_jobs << endl;
}

void print_cpu_info(vector<DirectedGraph::vertex_descriptor> run_on_proc) {
    cout << "=== CPU INFO ===" << endl;
    for (auto i = 0; i < num_proc; ++i) {
        cout << "CPU " << i << ": " << run_on_proc[i] << endl;
    }
    cout << "=== CPU INFO ===" << endl;
}

void hybrid_logic_blox_scheduler(size_t num_levels) {
    bool CPUfree = false;
    bool allCPUfree = false;
    vector<double> end_time_on_proc(num_proc, 0.0);
    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc, 0);

    if (scheduled_nodes.size() == 0) {
        cout << "No nodes scheduled." << endl;
        exit(1);
    }

    // Schedule first eight jobs on processors
    // TODO: fix a case when there is less than eight jobs initially available
    if (scheduled_nodes.size() < 8) {
        cout << "Less than eight initial nodes currently is not supported." << endl;
        exit(1);
    }

    // Scheduling first jobs on the CPUs
    for (auto i = 0; i < num_proc; ++i) {
        auto node = get_hybrid_scheduled_node();
        if (node == -1) {
            for (auto j = i; j < num_proc; ++j) {
                end_time_on_proc[j] = 0.0;
                run_on_proc[j] = sink_node;
            }
            break;
        }
        end_time_on_proc[i] = G[node].proc_time;
        run_on_proc[i] = node;
        if (scheduled_nodes.find(node) == scheduled_nodes.end() && node != -1) {
            cout << "ERROR: Trying to schedule node not from the scheduled set (1)." << endl;
            exit(1);
        }
        scheduled_nodes.erase(node);
        active_nodes_on_level[G[node].level].erase(node);
        running_nodes.insert(node);
    }

    while (scheduled_nodes.size() != 0 || !allCPUfree || !run_queue.empty()) {
        cout << "Keep scheduling..." << endl;
        if (!CPUfree) {
            cout << "Node is not free. Finishing one of the jobs..." << endl;
            finish_hybrid_earliest_job(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        }
        auto node = get_hybrid_scheduled_node();
        if (scheduled_nodes.find(node) == scheduled_nodes.end() && node != -1) {
            cout << "ERROR: Trying to schedule node not from the scheduled set (2)." << endl;
            exit(1);
        }
        while (node == -1) {
            cout << "Trying to schedule by finishing work..." << endl;
            finish_hybrid_earliest_job(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            cout << "Finished the job." << endl;
            node = get_hybrid_scheduled_node();
            if (scheduled_nodes.find(node) == scheduled_nodes.end() && node != -1) {
                cout << "ERROR: Trying to schedule node not from the scheduled set (3)." << endl;
                exit(1);
            }
            if (allCPUfree && (node == -1))
                break;
        }
        cout << "Breaked out of the loop." << endl;
        if (allCPUfree && (node == -1))
            continue;

        int cpu_num = -1;
        double min_time = numeric_limits<double>::max();
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] != sink_node)
                continue;
            if (end_time_on_proc[i] < min_time) {
                min_time = end_time_on_proc[i];
                cpu_num = i;
            }
        }
        if (cpu_num == -1) {
            cout << "ERROR: Scheduling error." << endl;
            exit(1);
        }
        run_on_proc[cpu_num] = node;
        end_time_on_proc[cpu_num] += G[node].proc_time;
        scheduled_nodes.erase(node);
        active_nodes_on_level[G[node].level].erase(node);
        running_nodes.insert(node);
        allCPUfree = false;
        CPUfree = false;
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] == sink_node) {
                CPUfree = true;
                break;
            }
        }
//        print_cpu_info(run_on_proc);
    }

    double max_time = -1.0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] > max_time) {
            max_time = end_time_on_proc[i];
        }
    }

    cout << "Total running time of schedule: " << max_time << endl;
    cout << "Number of executed jobs: " << num_executed_jobs << endl;
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

// TODO: we don't really need that function to be separate for the LB and SBU scheduler. Move it to the defs.hpp
// with some minor changes.
void lb_add_start_nodes() {
    ifstream input_file;
    input_file.open("startingNodes.dat");

    uint64_t node;
    string node_orig;

    auto num_active_nodes = 0;
    size_t num_lines_in_file = count_lines_in_file("startingNodes.dat");

    for (auto i = 0; i < num_lines_in_file; ++i) {
        input_file >> node_orig;
        G[vertices_mapping[node_orig]].is_active = true;
        num_active_nodes++;
        start_nodes.insert(vertices_mapping[node_orig]);
        if (G[vertices_mapping[node_orig]].is_predicate_node) {
            cout << "Starting node is a predicate node." << endl;
        } else {
            cout << "Starting node is NOT a predicate node." << endl;
        }
    }

#ifdef DEBUG_PRINT
    cout << "Number of starting nodes (from file): " << num_active_nodes << endl;
#endif

    input_file.close();
}

bool is_ancestor(DirectedGraph::vertex_descriptor child, DirectedGraph::vertex_descriptor parent) {
    lookup_node = child;
    lookup_node_found = false;
    ancestral_bfs vis;
    breadth_first_search(G, parent, visitor(vis));

    return lookup_node_found;
}

// TODO: currently this function is implemented in a very naive way so
// I may need to optimize it later using algorithm from LogicBlox
void calculate_ancestral_lists(vector<DirectedGraph::vertex_iterator> source_nodes) {
    ancestral_list.resize(num_vertices(G) + 1);

    // TODO: This functionality does not work yet (maybe this is wrong or the way the file is written
    if (parse_ancestor_list) {
        ifstream input_file;
        input_file.open("ancestorList.dat");

        uint64_t main_node;
        uint64_t node_to_insert;

        size_t num_lines_in_file = count_lines_in_file("ancestorList.dat");

        if (num_lines_in_file != num_vertices(G) + 1) {
            cout << "Number of nodes in the ancestor list file does not match number of nodes in the current graph." << endl;
            exit(EINVAL);
        }

        string str;
        while(getline(input_file, str)) {
            istringstream ss(str);
            ss >> main_node;
            while(ss >> node_to_insert) {
                ancestral_list[main_node].insert(node_to_insert);
            }
        }

        input_file.close();

        return ;
    }

    // First check that all source nodes have empty ancestral lists
    for (auto it : source_nodes) {
        if (ancestral_list[*it].size() != 0) {
            cout << "ERROR: ancestral list of source node " << *it << " is not empty." << endl;
            exit(1);
        }
    }

    add_vertex(G);
    auto root_node = num_vertices(G) - 1;
    for (auto it : source_nodes) {
        add_edge(root_node, *it, G);
    }

    // Run the loop twice as otherwise some nodes may be visited before BFS will visit their ancestors
    cout << "Running ancestral BFS..." << endl;
    auto iter_count = 0;
    while (!ancestor_list_converged) {
        ancestor_list_converged = true;
        merge_ancestral_info_bfs vis;
        breadth_first_search(G, root_node, visitor(vis));
        iter_count++;
    }
    cout << "Complete. It took " << iter_count << " iterations." << endl;

    clear_out_edges(root_node, G);

    ofstream output_file;
    output_file.open("ancestorList.dat");

    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        output_file << *vit;
        for (auto node : ancestral_list[*vit]) {
            output_file << " " << node;
        }
        output_file << endl;
    }

    output_file.close();
}
///////////////////////////////////////////////////////////////////////
void add_to_busy_set(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(node, G);
    for (; eit != eend; ++eit) {
        if (!G[target(*eit, G)].is_predicate_node) {
            cout << "ERROR: Trying to add non-predicate node to the busy set." << endl;
            exit(1);
        }
        busy_set.insert(target(*eit, G));
    }
}

void finish_earliest_job_helper(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
    int cpu_num = -1;
    double min_time = numeric_limits<double>::max();
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            continue;
        if (end_time_on_proc[i] < min_time) {
            min_time = end_time_on_proc[i];
            cpu_num = i;
        }
    }

    if (cpu_num == -1) {
        cout << "ERROR: No job to finish." << endl;
        exit(1);
    }

    running_nodes.erase(run_on_proc[cpu_num]);

    // Once the active job finishes, the end time on all idle cpus
    // also should be updated
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            end_time_on_proc[i] = min_time;
    }

    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(run_on_proc[cpu_num], G);
    for (; eit != eend; ++eit) {
        if (!G[target(*eit, G)].is_predicate_node) {
            cout << "ERROR: Tryong to delete non-predicate node from the busy set." << endl;
            exit(1);
        }
        if (busy_set.find(target(*eit, G)) == busy_set.end()) {
            cout << "ERROR: Busy set does not contain node to be removed." << endl;
            exit(1);
        }
        busy_set.erase(target(*eit, G));
        for (auto el : busy_watch_set[target(*eit, G)]) {
            G[el].in_busy_watch = false;
            activation_queue.insert(el);
        }
        busy_watch_set[target(*eit, G)].clear();

        DirectedGraph::out_edge_iterator eit2, eend2;
        tie(eit2, eend2) = out_edges(target(*eit, G), G);
        for (; eit2 != eend2; ++eit2) {
            if (G[target(*eit2, G)].is_active && !G[target(*eit2, G)].in_busy_watch) {
                activation_queue.insert(target(*eit2, G));
            }
        }
    }

    // TODO: potentially add later to account for scheduling overhead
//    current_timestamp = end_time_on_proc[cpu_num];

    G[run_on_proc[cpu_num]].is_executed = true;
    if (!G[run_on_proc[cpu_num]].is_active) {
        cout << "Executed non-active job." << endl;
        exit(1);
    }

    run_on_proc[cpu_num] = sink_node;
    CPUfree = true;
    num_executed_jobs++;

    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] != sink_node)
            return ;
    }

    allCPUfree = true;
}

DirectedGraph::vertex_descriptor get_node_helper(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
    sched_start = chrono::high_resolution_clock::now();
    DirectedGraph::vertex_descriptor node = -1;

    if (activation_queue.empty()) {
//        sched_end = chrono::high_resolution_clock::now();
//        scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

        return -1;
    }

    auto el = activation_queue.begin();
    bool should_be_moved = false;
    while (el != activation_queue.end()) {
        should_be_moved = false;
        for (auto busy_predicate_node : busy_set) {
            if (ancestral_list[*el].find(busy_predicate_node) != ancestral_list[*el].end()) {
                busy_watch_set[busy_predicate_node].insert(*el);
                should_be_moved = true;
                break;
            }
        }
        if (should_be_moved) {
            set<DirectedGraph::vertex_descriptor, decltype(cmp)>::iterator it = el++;
            G[*it].in_busy_watch = true;
            activation_queue.erase(*it);
        } else {
            node = *el;
            break;
        }
    }

    sched_end = chrono::high_resolution_clock::now();
    scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

    return node;
}

void lb_scheduler_helper() {
    bool CPUfree = true;
    bool allCPUfree = false;
    vector<double> end_time_on_proc(num_proc, 0.0);
    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc, 0);

//    cout << "Activation queue size (init): " << activation_queue.size() << endl;
    for (auto node : start_nodes) {
        if (!G[node].is_predicate_node) {
            cout << "ERROR: Starting node is NOT a predicate node." << endl;
            exit(1);
        }
        lb_activate_children_helper(node);
    }
//    cout << "Activation queue size (after init): " << activation_queue.size() << endl;


    if (activation_queue.size() == 0) {
        cout << "No nodes scheduled." << endl;
        exit(1);
    }

    // Schedule first eight jobs on processors
    // TODO: fix a case when there is less than eight jobs initially available
    if (activation_queue.size() < 8) {
        cout << "Less than eight initial nodes currently is not supported." << endl;
        exit(1);
    }

    // Init CPUs
    // TODO: this can be moved to the initialization of these variables
    for (auto i = 0; i < num_proc; ++i) {
        end_time_on_proc[i] = 0.0;
        run_on_proc[i] = sink_node;
    }

    while (activation_queue.size() != 0 || !allCPUfree) {
        if (!CPUfree) {
            finish_earliest_job_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        }
        auto node = get_node_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        while (node == -1) {
            finish_earliest_job_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            node = get_node_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            if (allCPUfree && (node == -1))
                break;
        }
        if (allCPUfree && (node == -1))
            continue;

        int cpu_num = -1;
        double min_time = numeric_limits<double>::max();
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] != sink_node)
                continue;
            if (end_time_on_proc[i] < min_time) {
                min_time = end_time_on_proc[i];
                cpu_num = i;
            }
        }
        run_on_proc[cpu_num] = node;
        end_time_on_proc[cpu_num] += G[node].proc_time;
        activation_queue.erase(node);
        running_nodes.insert(node);
        add_to_busy_set(node);
        allCPUfree = false;
        CPUfree = false;
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] == sink_node) {
                CPUfree = true;
                break;
            }
        }
    }

    double max_time = -1.0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] > max_time) {
            max_time = end_time_on_proc[i];
        }
    }

    cout << "Total running time of schedule: " << max_time << endl;
//    cout << "Number of executed jobs: " << num_executed_jobs << endl;
}
///////////////////////////////////=========================== HYBRID
void add_to_busy_set_hybrid(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(node, G);
    for (; eit != eend; ++eit) {
        if (!G[target(*eit, G)].is_predicate_node) {
            cout << "ERROR: Trying to add non-predicate node to the busy set." << endl;
            exit(1);
        }
        busy_set.insert(target(*eit, G));
    }
}

void finish_earliest_job_hybrid_helper(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
    int cpu_num = -1;
    double min_time = numeric_limits<double>::max();
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            continue;
        if (end_time_on_proc[i] < min_time) {
            min_time = end_time_on_proc[i];
            cpu_num = i;
        }
    }

    if (cpu_num == -1) {
        cout << "ERROR: No job to finish." << endl;
        exit(1);
    }

    running_nodes.erase(run_on_proc[cpu_num]);

    // Once the active job finishes, the end time on all idle cpus
    // also should be updated
    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] == sink_node)
            end_time_on_proc[i] = min_time;
    }

    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(run_on_proc[cpu_num], G);
    for (; eit != eend; ++eit) {
        if (!G[target(*eit, G)].is_predicate_node) {
            cout << "ERROR: Tryong to delete non-predicate node from the busy set." << endl;
            exit(1);
        }
        if (busy_set.find(target(*eit, G)) == busy_set.end()) {
            cout << "ERROR: Busy set does not contain node to be removed." << endl;
            exit(1);
        }
        busy_set.erase(target(*eit, G));
        for (auto el : busy_watch_set[target(*eit, G)]) {
            G[el].in_busy_watch = false;
            activation_queue.insert(el);
        }
        busy_watch_set[target(*eit, G)].clear();

        DirectedGraph::out_edge_iterator eit2, eend2;
        tie(eit2, eend2) = out_edges(target(*eit, G), G);
        for (; eit2 != eend2; ++eit2) {
            if (G[target(*eit2, G)].is_active && !G[target(*eit2, G)].in_busy_watch) {
                activation_queue.insert(target(*eit2, G));
            }
            if (G[target(*eit2, G)].is_active) {
                active_nodes_on_level[G[target(*eit, G)].level].insert(target(*eit2, G));
            }
        }
    }

    if (G[run_on_proc[cpu_num]].level == current_level)
        num_run_nodes_on_level--;

    if (run_queue.empty() && num_run_nodes_on_level == 0) {
        if (current_level >= find_depth()) {
            ;
        } else {
            while(num_run_nodes_on_level == 0) {
                current_level++;
                cout << "This is level: " << current_level << endl;
                if (current_level >= find_depth())
                    break;
                for (auto el : active_nodes_on_level[current_level]) {
                    run_queue.insert(el);
                }
                num_run_nodes_on_level = run_queue.size();
                cout << "Moving levels. New run queue size: " << run_queue.size() << endl;
                cout << "This is level: " << current_level << endl;
            }
        }
    }
    if (num_run_nodes_on_level < 0) {
        cout << "negative level" << endl;
        exit(1);
    }
    active_nodes_on_level[G[run_on_proc[cpu_num]].level].erase(run_on_proc[cpu_num]);

    // TODO: potentially add later to account for scheduling overhead
//    current_timestamp = end_time_on_proc[cpu_num];

    G[run_on_proc[cpu_num]].is_executed = true;
    if (!G[run_on_proc[cpu_num]].is_active) {
        cout << "Executed non-active job." << endl;
        exit(1);
    }

    run_on_proc[cpu_num] = sink_node;
    CPUfree = true;
    num_executed_jobs++;

    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] != sink_node)
            return ;
    }

    allCPUfree = true;
}

DirectedGraph::vertex_descriptor get_node_hybrid_helper(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
    sched_start = chrono::high_resolution_clock::now();
    DirectedGraph::vertex_descriptor node = -1;

    if (!run_queue.empty()) {
        node = *(run_queue.begin());
        sched_end = chrono::high_resolution_clock::now();
        scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

        cout << "Node is coming from levels..." << endl;
        return node;
    }

    if (activation_queue.empty()) {
//        sched_end = chrono::high_resolution_clock::now();
//        scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

        return -1;
    }

    auto el = activation_queue.begin();
    bool should_be_moved = false;
    while (el != activation_queue.end()) {
        should_be_moved = false;
        for (auto busy_predicate_node : busy_set) {
            if (ancestral_list[*el].find(busy_predicate_node) != ancestral_list[*el].end()) {
                busy_watch_set[busy_predicate_node].insert(*el);
                should_be_moved = true;
                break;
            }
        }
        if (should_be_moved) {
            set<DirectedGraph::vertex_descriptor, decltype(cmp)>::iterator it = el++;
            G[*it].in_busy_watch = true;
            activation_queue.erase(*it);
        } else {
            node = *el;
            break;
        }
    }

    sched_end = chrono::high_resolution_clock::now();
    scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

    return node;
}

void lb_hybrid_activate_children_helper(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(node, G);
    // TODO: add a check against busy watch set
    for (; eit != eend; ++eit) {
        if (G[target(*eit, G)].is_active) {
            if (G[target(*eit, G)].is_predicate_node) {
                cout << "ERROR: Trying to activate predicate node." << endl;
                exit(1);
            }
            activation_queue.insert(target(*eit, G));
            active_nodes_on_level[G[target(*eit, G)].level].insert(target(*eit, G));
        }
    }

    for (auto el : active_nodes_on_level[1]) {
        run_queue.insert(el);
    }
    num_run_nodes_on_level = run_queue.size();
}

void lb_hybrid_scheduler_helper() {
    bool CPUfree = true;
    bool allCPUfree = false;
    vector<double> end_time_on_proc(num_proc, 0.0);
    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc, 0);

    for (auto node : start_nodes) {
        if (!G[node].is_predicate_node) {
            cout << "ERROR: Starting node is NOT a predicate node." << endl;
            exit(1);
        }
        lb_hybrid_activate_children_helper(node);
    }


    if (activation_queue.size() == 0) {
        cout << "No nodes scheduled." << endl;
        exit(1);
    }

    // Schedule first eight jobs on processors
    // TODO: fix a case when there is less than eight jobs initially available
    if (activation_queue.size() < 8) {
        cout << "Less than eight initial nodes currently is not supported." << endl;
        exit(1);
    }

    // Init CPUs
    // TODO: this can be moved to the initialization of these variables
    for (auto i = 0; i < num_proc; ++i) {
        end_time_on_proc[i] = 0.0;
        run_on_proc[i] = sink_node;
    }

    while (activation_queue.size() != 0 || !allCPUfree) {
        if (!CPUfree) {
            finish_earliest_job_hybrid_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        }
        auto node = get_node_hybrid_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        while (node == -1) {
            finish_earliest_job_hybrid_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            node = get_node_hybrid_helper(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            if (allCPUfree && (node == -1))
                break;
        }
        if (allCPUfree && (node == -1))
            continue;

        int cpu_num = -1;
        double min_time = numeric_limits<double>::max();
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] != sink_node)
                continue;
            if (end_time_on_proc[i] < min_time) {
                min_time = end_time_on_proc[i];
                cpu_num = i;
            }
        }
        run_on_proc[cpu_num] = node;
        end_time_on_proc[cpu_num] += G[node].proc_time;
        activation_queue.erase(node);
        run_queue.erase(node);
        running_nodes.insert(node);
        add_to_busy_set_hybrid(node);
        allCPUfree = false;
        CPUfree = false;
        for (auto i = 0; i < num_proc; ++i) {
            if (run_on_proc[i] == sink_node) {
                CPUfree = true;
                break;
            }
        }
    }

    double max_time = -1.0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] > max_time) {
            max_time = end_time_on_proc[i];
        }
    }
}
///////////////////////////////////===========================
void lb_scheduler(bool trace_enabled) {
    cout << "Starting LogicBlox scheduler..." << endl;
    cout << "vertices: " << num_vertices(G) << endl;
    cout << "edges: " << num_edges(G) << endl;

    vector<DirectedGraph::vertex_iterator> source_nodes;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);

    uint64_t many_count = 0;
    uint64_t one_count = 0;
    uint64_t zero_count = 0;

    busy_watch_set.resize(num_vertices(G) + 1);

    lb_init_graph(source_nodes);
//    rebuild_graph();
//    find_levels_ds(G_tmp);

    find_levels_ds(G);

    // Adding levels here for the hybrid scheduler
    size_t num_levels = find_depth();
    active_nodes_on_level.resize(num_levels + 1);

    cout << "Number of levels: " << num_levels << endl;
    exit(1);

    // TODO: add extensive graph/trace check function so to remove them from the various functions

    if (trace_enabled) {
        lb_add_start_nodes();
        assign_jobs_from_trace();

        cout << "Traces are enabled." << endl;
        cout << "vertices: " << num_vertices(G) << endl;
        cout << "edges: " << num_edges(G) << endl;
        cout << "vertices: " << num_vertices(G_tmp) << endl;
        cout << "edges: " << num_edges(G_tmp) << endl;
        cout << "Number of source nodes: " << source_nodes.size() << endl;
    } else {
        cout << "LogicBlox scheduler without traces is not supported yet." << endl;
        exit(EINVAL);
    }

    /*
    uint64_t num_start_nodes = 0;
    for (auto node : source_nodes) {
        if (G[*node].is_active) {
            num_start_nodes--;
            lb_activate_children(*node);
        }
    }

    if (num_start_nodes != 0) {
        cout << "One of the starting nodes was not source nodes." << endl;
        exit(1);
    }
    */

    /*
    uint64_t num_act = 0;
    for (; vit != vend; ++vit) {
        if (G[*vit].is_active) num_act++;
    }
    cout << "ACTIVE JOBS: " << num_act << endl;
    */

    cout << "Calculating ancestral lists..." << endl;
    calculate_ancestral_lists(source_nodes);
    cout << "Done calculating ancestral lists." << endl;
    add_sink_node();

    cout << "Size of the start set: " << start_nodes.size() << endl;

//    auto sm_sched_start = chrono::high_resolution_clock::now();
    auto start_t = omp_get_wtime();
    for (auto ij = 0; ij < 1; ij++)
        lb_scheduler_helper();
//        lb_hybrid_scheduler_helper();
    auto end_t = omp_get_wtime();
//    auto sm_sched_end = chrono::high_resolution_clock::now();
//    auto sm_scheduling_overhead = chrono::duration_cast<chrono::microseconds>(sm_sched_end - sm_sched_start);

//    cout << "Total time: " << sm_scheduling_overhead.count() / 10.0<< endl;
    cout << "Total time: " << (double)(end_t - start_t) / 1.0 << endl;
    cout << "Number of executed jobs: " << num_executed_jobs / 1 << endl;

    cout << "Total scheduling overhead: " << scheduling_overhead.count() << endl;
    cout << "Total scheduling overhead: " << scheduling_overhead.count() / 1000 << endl;

    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        if (G[*vit].is_active && !G[*vit].is_executed && !G[*vit].is_predicate_node) {
            cout << "NOT EXECUTED: " << *vit << endl;
            cout << "LEVEL: " << G[*vit].level << endl;
        }
    }

    return ;

    cout << "Active jobs on level 0: " << active_nodes_on_level[0].size() << endl;
    cout << "Active jobs on level 1: " << active_nodes_on_level[1].size() << endl;

    num_run_nodes_on_level = active_nodes_on_level[1].size();
    current_level = 1;
    for (auto node : active_nodes_on_level[1]) {
        run_queue.insert(node);
        scheduled_nodes.erase(node);
    }


    current_timestamp = 0.0;
    cout << "Initial scheduling overhead (should be zero): " << scheduling_overhead.count() << endl;
    for (; vit != vend; ++vit)
        if (G[*vit].seq_number == 5001)
            cout << "FOUDN:: " << G[*vit].level << endl;
    cout << "FOUDN:: " << G[1282].level << endl;
    exit(1);

    logic_blox_scheduler();

    if (num_levels != find_depth()) {
        cout << "ERROR: Number of levels changed." << endl;
        exit(1);
    }
//    hybrid_logic_blox_scheduler(num_levels + 1);

    if (print_executed_nodes)
        executed_nodes_file.close();

    cout << "Total scheduling overhead: " << scheduling_overhead.count() << endl;
    cout << "Total scheduling overhead: " << scheduling_overhead.count() / 1000 << endl;

    cout << "Done." << endl;
}
