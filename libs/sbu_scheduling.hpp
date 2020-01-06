vector<set<DirectedGraph::vertex_descriptor> > running_nodes_on_level;
vector<DirectedGraph::vertex_descriptor> running_vec;
auto sbu_total_t = omp_get_wtime();
uint64_t bfs_active_nodes;
set<DirectedGraph::vertex_descriptor> bfs_active_set;

set<DirectedGraph::vertex_descriptor> subgraph_nodes;
set<DirectedGraph::edge_descriptor> subgraph_edges;
auto vertices_count = 0;
auto edges_count = 0;

struct custom_bfs_visitor : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        if (G[target(e, G)].level <= G[source(e, G)].level) {
            if (G[target(e, G)].is_predicate_node) {
                // Uncomment one below and comment another to prevent 
                // having separate level for predicate nodes
//                G[target(e, G)].level = G[source(e, G)].level;
                G[target(e, G)].level = G[source(e, G)].level + 1;
            } else {
                G[target(e, G)].level = G[source(e, G)].level + 1;
            }
        }
    }
};

DirectedGraph::vertex_descriptor node_in_search;
bool node_in_search_found;

struct searching_bfs : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        if (target(e, G) == node_in_search) {
            node_in_search_found = true;
        }
    }
};

struct subgraph_bfs : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        subgraph_nodes.insert(target(e, g));
        subgraph_edges.insert(e);
        edges_count++;
    }
};

struct counting_bfs : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        if (G[target(e, G)].is_active) {
            bfs_active_nodes++;
            bfs_active_set.insert(target(e, G));
        }
    }
};

void sbu_init_graph(vector<DirectedGraph::vertex_iterator> &source_nodes) {
#ifdef DEBUG_PRINT
    cout << "Initializing vertex properties and constructing set of source nodes..." << endl;
#endif
    uint64_t source_unit_nodes = 0;
    num_executed_jobs = 0;

    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    DirectedGraph::in_edge_iterator inedgeIt, inedgeEnd;
    tie(vertexIt, vertexEnd) = vertices(G);

    for (; vertexIt != vertexEnd; ++vertexIt) {
        G[*vertexIt].level = -1;
        G[*vertexIt].level_mod = -1;
        G[*vertexIt].is_active = false;
        G[*vertexIt].is_scheduled = false;
        G[*vertexIt].proc_time = 0;
        G[*vertexIt].in_degree = in_degree(*vertexIt, G);
        G[*vertexIt].prunned = false;
        if (G[*vertexIt].in_degree == 0) {
            source_nodes.push_back(vertexIt);
            G[*vertexIt].is_source_node = true;
            if (!G[*vertexIt].is_predicate_node) {
                cout << "WARNING: Source unit node.\r";
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

void populate_levels() {
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++ vertexIt) {
        if (!G[*vertexIt].prunned) {
            levels[G[*vertexIt].level].push_back(vertexIt);
        }
    }
}

void check_dog() {
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++ vertexIt) {
        if (G[*vertexIt].level == 0 && G[*vertexIt].is_scheduled) {
            cout << "FARRKY" << endl;
            exit(1);
        }
    }
}

void run_sbu_job(int level) {
    this_thread::sleep_for(chrono::milliseconds(10));
    lock_guard<mutex> lock(jobs_mutex);
    auto node = running_nodes_on_level[level].begin();
    DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
    tie(edgeIt, edgeEnd) = out_edges(*node, G);
    for (; edgeIt != edgeEnd; ++edgeIt) {
        if (G[target(*edgeIt, G)].is_active) {
            active_nodes_on_level[G[target(*edgeIt, G)].level].insert(target(*edgeIt, G));
        }
    }
    running_nodes_on_level[level].erase(*node);
}

void find_levels_ds(DirectedGraph graph) {
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(graph);
    for (; vit != vend; ++vit) {
        if (graph[*vit].level_mod != -1) {
            cout << "ERROR: wrong initial level value." << endl;
            exit(1);
        }
        graph[*vit].in_degree = in_degree(*vit, graph);
    }
    uint64_t max_in_degree = 0;
    tie(vit, vend) = vertices(graph);
    for (; vit != vend; ++vit) {
        if (graph[*vit].in_degree > max_in_degree) {
            max_in_degree = graph[*vit].in_degree;
        }
    }
    vector<set<DirectedGraph::vertex_descriptor> > leveling_ds(max_in_degree + 1);

    tie(vit, vend) = vertices(graph);
    for (; vit != vend; ++vit) {
        if (graph[*vit].is_source_node) {
            leveling_ds[0].insert(*vit);
            if (graph[*vit].in_degree != 0) {
                cout << "ERROR: in-degree of source node is not zero." << endl;
                exit(1);
            }
        } else {
            leveling_ds[graph[*vit].in_degree].insert(*vit);
            if (graph[*vit].in_degree == 0 && !graph[*vit].is_predicate_node) {
                cout << "ERROR: unit node with in-degree zero." << endl;
            }
        }
    }
    size_t current_level = 0;
    size_t ds_size = 1;
    while (ds_size != 0) {
        ds_size = 0;
        for (auto node = leveling_ds[0].begin(); node != leveling_ds[0].end(); ) {
            graph[*node].level_mod = current_level;
            clear_out_edges(*node, graph);
            node = leveling_ds[0].erase(node);
        }
        for (auto i = 1; i < max_in_degree + 1; ++i) {
            for (auto node = leveling_ds[i].begin(); node != leveling_ds[i].end(); ) {
                graph[*node].in_degree = in_degree(*node, graph);
                if (graph[*node].in_degree != i) {
                    leveling_ds[graph[*node].in_degree].insert(*node);
                    node = leveling_ds[i].erase(node);
                } else {
                    node++;
                }
            }
        }
        for (auto i = 0; i < max_in_degree + 1; ++i) {
            ds_size += leveling_ds[i].size();
        }
        if (ds_size != 0 && leveling_ds[0].size() == 0) {
            cout << "ERROR: no nodes with in-degree zero in non-empty graph." << ds_size << endl;
            exit(1);
        }
        current_level++;
    }

    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        if (graph[*vit].level_mod < 0) {
            cout << "ERROR: node without assigned level." << endl;
            exit(1);
        }
        G[*vit].level = graph[*vit].level_mod;
    }

    cout << "Max of level: " << current_level - 1 << endl;
}

void run_sbu_fake_job(int level) {
//    lock_guard<mutex> lock(jobs_mutex);
    auto node = running_nodes_on_level[level].begin();
    DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
    tie(edgeIt, edgeEnd) = out_edges(*node, G);
    for (; edgeIt != edgeEnd; ++edgeIt) {
        if (G[target(*edgeIt, G)].is_active) {
            active_nodes_on_level[G[target(*edgeIt, G)].level].insert(target(*edgeIt, G));
        }
    }
    running_nodes_on_level[level].erase(*node);
}

void schedule_sbu_job(set<uint64_t> &active_jobs, int level) {
//    lock_guard<mutex> lock(jobs_mutex);
    auto node = active_jobs.begin();
    auto start_t = omp_get_wtime();
    running_nodes_on_level[level].insert(*node);
    active_jobs.erase(*node);
    auto end_t = omp_get_wtime();
    sbu_total_t += (end_t - start_t);
    run_sbu_fake_job(level);
    // TODO: seems like lock issues. it should be unlocked before calling run_job fxn or
    // can go into deadlock. fix it.
//    thread thr = thread(&run_job, level);
//    thr.detach();
}

double calculate_time(vector<vector<DirectedGraph::vertex_iterator> > &levels, uint64_t num_levels) {
    double result = 0;
    uint64_t counter = 0;
    for (auto i = 1; i < num_levels; ++i) {
        double num_active_nodes = 0;
        for (auto j : levels[i]) {
            if (G[*j].is_active) {
                num_active_nodes++;
            }
        }
        result += ceil(num_active_nodes / num_proc);
        counter += num_active_nodes;
    }
#ifdef DEBUG_PRINT
    cout << "Total number of jobs: " << counter << endl;
#endif
    return result;
}

DirectedGraph::vertex_descriptor next_node;
DirectedGraph::vertex_descriptor guilty_node;
bool next_node_found;

struct lookahead_bfs_visitor : default_bfs_visitor {
    void examine_edge(const DirectedGraph::edge_descriptor &e, const DirectedGraph &g) const {
        if (G[target(e, G)].level > G[next_node].level) {
            return ;
        }
        if (target(e, G) == next_node) {
            next_node_found = true;
            guilty_node = source(e, G);
        }
    }
};

DirectedGraph::vertex_descriptor lookahead(uint64_t level, uint64_t num_levels) {
    lookahead_bfs_visitor vis;
    for (auto i = level + 1; i < std::min(num_levels, level + 15); ++i) {
        if (active_nodes_on_level[i].empty()) {
            continue;
        }
        for (auto elem : active_nodes_on_level[i]) {
            if (!G[elem].is_scheduled) {
                continue;
            }
            next_node_found = false;
            next_node = elem;
            for (auto el : running_vec) {
                breadth_first_search(G, el, visitor(vis));
            
            }
            if (!next_node_found) {
                active_nodes_on_level[i].erase(next_node);
                return next_node;
            }
        }
    }

    return 0;
}

void sbu_activate_children(DirectedGraph::vertex_descriptor node) {
    DirectedGraph::out_edge_iterator eit, eend;
    tie(eit, eend) = out_edges(node, G);
    for (; eit != eend; ++eit) {
        if (G[target(*eit, G)].is_active) {
            G[target(*eit, G)].is_scheduled = true;
            if (G[target(*eit, G)].level == 0) {
                cout << "CRITICAL ERROR" << endl;
                exit(1);
            }
            scheduled_nodes.insert(target(*eit, G));
        }
    }
}

// This experiment simulates level-based algorithm with internal parallelism. It is done by
// running tasks one-by-one using all available processors
void fully_parallel_experiment(uint64_t num_levels) {
#ifdef DEBUG_PRINT
    cout << "Starting SBU scheduler." << endl;
#endif
    double processing_time = 0;
    auto total_jobs = 0;

    for (auto level = 1; level < num_levels; ++level) {
        cout << "Processing level " << level << " with " << active_nodes_on_level[level].size() << " active nodes." << endl;
        if (active_nodes_on_level[level].size() == 0) {
            cout << "Level is empty. Skipping..." << endl;
            continue;
        }
        for (auto node : active_nodes_on_level[level]) {
            processing_time += 1.0 * G[node].proc_time / num_proc;
            total_jobs++;
        }
    }

    cout << "===== Fully parallel experiment =====" << endl;
    cout << "Number of tasks processed: " << total_jobs << endl;
    cout << "Processing time: " << processing_time << endl;

    return ;
}

void run_sbu_scheduler(uint64_t num_levels) {
#ifdef DEBUG_PRINT
    cout << "Starting SBU scheduler." << endl;
#endif
    sbu_total_t = 0;
    double processing_time = 0;
    auto total_jobs = 0;
    int num_running = 0;
    uint64_t cur_level;

    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc);
    vector<double> end_time_on_proc(num_proc);
    for (auto i = 0; i < num_proc; ++i) {
        run_on_proc[i] = sink_node;
        end_time_on_proc[i] = 0;
    }

    uint64_t init_num_nodes = (active_nodes_on_level[1].size() < num_proc) ? active_nodes_on_level[1].size() : num_proc;
    if (init_num_nodes == 0) {
        cout << "WARNING: level 1 is empty." << endl;
        sleep(3);
    }
    for (auto i = 0; i < init_num_nodes; ++i) {
        auto node = active_nodes_on_level[1].begin();
        run_on_proc[i] = *node;
        scheduled_nodes.insert(*node);
        end_time_on_proc[i] = G[*node].proc_time;
        active_nodes_on_level[1].erase(*node);
    }

    for (auto level = 1; level < num_levels; ++level) {
        cout << "Processing level " << level << " with " << active_nodes_on_level[level].size() << " active nodes." << endl;
        if (active_nodes_on_level[level].size() == 0) {
            cout << "Level is empty. Skipping..." << endl;
            continue;
        }
        for (auto node : active_nodes_on_level[level]) {
            if (G[node].is_predicate_node) {
                cout << "predicate node" << endl;
                exit(1);
            }
            if (!G[node].is_active) {
                cout << "not active" << endl;
                exit(1);
            }
            if (node == 61124) {
                cout << "is it active? " << G[node].is_active << endl;
                cout << "in degree? " << in_degree(node, G) << endl;
            }

            /* We should not check for that as parents nodes on the previous level
             * may still be running and thus have not activated any node on the next level
            if (!G[node].is_scheduled) {
                cout << "ERROR: found inactive job on current level " << level << " with ID: " << node << endl;

                DirectedGraph::in_edge_iterator it, ite;
                tie(it, ite) = in_edges(node, G);
                for ( ; it != ite; ++it) {
                    if (G[source(*it, G)].is_active) {
                        cout << "scheduled? " << G[source(*it, G)].is_scheduled << endl;
                        cout << "id: " << source(*it, G) << endl;
                        cout << "level: " << G[source(*it, G)].level << endl;
                    }
                }

                exit(1);
            }
            */
        }
        while (!active_nodes_on_level[level].empty()) {
            cur_level = level;
            auto node = active_nodes_on_level[level].begin();

            /* It is a heuristic that trying to schedule node with max out-degree first*/
            auto max_out_degree_node = node;
            //uint64_t max_out_degree = out_degree(*node, G);
            for (auto it = node; it != active_nodes_on_level[level].end(); ++it) {
                if (out_degree(*it, G) > out_degree(*max_out_degree_node, G)) {
                    max_out_degree_node = it;
                }
            }
            node = max_out_degree_node;
            /**/


            if (level > 1) {
                while (1) {
                    // No need to check if the node is scheduled as it may depend on yet running nodes
                    // meaning being not scheduled yet. If it depends on the running node it will fail
                    // safety check even though it will be a candidate to run. Situation when it is checked
                    // and set as a safe to run should not occur due to the leveling nature and initial checks
                    // during the setup. If number of nodes scheduled to run will be correct than we assume
                    // that this part is indeed correct.
                    node = active_nodes_on_level[level].begin();

                    /* It is a heuristic that trying to schedule node with max out-degree first*/
                    max_out_degree_node = node;
                    //uint64_t max_out_degree = out_degree(*node, G);
                    for (auto it = node; it != active_nodes_on_level[level].end(); ++it) {
                        if (out_degree(*it, G) > out_degree(*max_out_degree_node, G)) {
                            max_out_degree_node = it;
                        }
                    }
                    node = max_out_degree_node;
                    /**/

                    double max = 0;
                    auto max_index = 0;
                    bool exhaust = true;
                    // Finding processor that running job with latest finishing time
                    for (auto i = 0; i < num_proc; ++i) {
                        if (end_time_on_proc[i] > max) {
                            max = end_time_on_proc[i];
                            max_index = i;
                        }
                    }
                    double min = max;
                    auto min_index = max_index;
                    // Finding processor that running job with earliest finishing time
                    // but excluding processors that weren't assigned any job yet
                    // TODO: why exclude 0?
                    for (auto i = 0; i < num_proc; ++i) {
                        if (end_time_on_proc[i] != 0 && end_time_on_proc[i] < min) {
                            min = end_time_on_proc[i];
                            min_index = i;
                        }
                    }

                    // If there is only one active node on the level then we check
                    // if it is safe to run (against all scheduled nodes) and if yes,
                    // then go straight to the scheduling it
                    if (active_nodes_on_level[level].size() == 1) {
#ifdef DBG_PRINT
                        cout << "Only one active node." << endl;
#endif
                        lookahead_bfs_visitor vis;
                        next_node_found = false;
                        next_node = *node;
//                        for (auto el : run_on_proc) {
//                        TODO:
//                        - instead of doing BFS check interval lists of scheduled nodes
                        for (auto el : scheduled_nodes) {
                            if (el == *node) {
                                continue;
                            }
                            breadth_first_search(G, el, visitor(vis));
                        }
                        if (!next_node_found) {
#ifdef DBG_PRINT
                            cout << "Node " << *node << " is safe to run (single node on level)." << endl;
#endif
                            exhaust = false;
                            break;
                        }
                    }

                    cur_level = level;
                    for (auto level2 = cur_level + 1; level2 < std::min(num_levels, cur_level + lookahead_levels); ++level2) {
                        for (node = active_nodes_on_level[level2].begin(); node != active_nodes_on_level[level2].end(); ++node) {
#ifdef DBG_PRINT
                            cout << "Performing look-ahead. Current depth: " << level2 - cur_level << endl;
#endif
                            // We check here if the node is scheduled as now we are trying to find one scheduled to run.
                            // If it is on the other level it definitely should be scheduled to run.
                            // TODO: implement heuristic of choosing node with higher out degree
                            // probably mark node as visited. If this works than better to use priority queue as
                            // a data structure
                            if (!G[*node].is_scheduled) {
                                continue;
                            }
                            lookahead_bfs_visitor vis;
                            next_node_found = false;
                            next_node = *node;
                            guilty_node = -1;
    //                        for (auto el : run_on_proc) {
    //                        TODO:
    //                        - instead of doing BFS check interval lists of scheduled nodes
                            for (auto el : scheduled_nodes) {
                                if (el == *node) {
                                    continue;
                                }
                                breadth_first_search(G, el, visitor(vis));
                                if (guilty_node != -1) {
                                    break;
                                }
                            }
                            if (!next_node_found) {
#ifdef DBG_PRINT
                                cout << "Node " << *node << " is safe to run." << endl;
#endif
                                exhaust = false;
                                cur_level = level2;
                                break;
                            }
                        }
                        if (!exhaust) {
                            break;
                        }
                    }
                    if (exhaust) {
#ifdef DBG_PRINT
                        cout << "===== Current running nodes =====" << endl;
                        for (auto j = 0; j < num_proc; ++j) {
                            cout << "Node ID:\t" << run_on_proc[j] << ".\tEnd time: " << end_time_on_proc[j] << endl;
                        }
                        cout << "Not found job to schedule. Finishing earliest one." << endl;
                        cout << "Current level is: " << level << endl;
                        cout << "Number of scheduled nodes: " << scheduled_nodes.size() << endl;
#endif        
                        double min2 = max;
                        auto min_index2 = max_index;
                        bool all_sink_nodes = true;
                        for (auto i = 0; i < num_proc; ++i) {
                            if (end_time_on_proc[i] != 0 && end_time_on_proc[i] <= min2 && run_on_proc[i] != num_vertices(G) - 1) {
                                min2 = end_time_on_proc[i];
                                min_index2 = i;
                                all_sink_nodes = false;
                            }
                        }

                        if (all_sink_nodes) {
                            node = active_nodes_on_level[level].begin();
#ifdef DBG_PRINT
                            cout << "Number of active jobs: " << active_nodes_on_level[level].size() << endl;
                            cout << "Last node attempted to be scheduled: " << *node << endl;
#endif
                            lookahead_bfs_visitor vis;
                            next_node_found = false;
                            next_node = *node;
                            guilty_node = -1;
                            auto el_tmp = guilty_node;
                            for (auto el : scheduled_nodes) {
                                if (el == *node) {
                                    continue;
                                }
                                breadth_first_search(G, el, visitor(vis));
                                if (guilty_node != -1) {

                                    el_tmp = el;
                                    cout << "Are they same? " << (el_tmp == guilty_node) << endl;
                                    cout << "EL Temporary: " << el_tmp << endl;
                                    break;
                                }
                            }
                            if (next_node_found) {
                                cout << "Guilty node: " << guilty_node << endl;
                                cout << "Our node level: " << G[*node].level << endl;
                                cout << "Guilty node level: " << G[guilty_node].level << endl;
                                cout << "is it active and scheduled? " << G[guilty_node].is_active << " " << G[guilty_node].is_scheduled << endl;
                                cout << "contains? " << (scheduled_nodes.find(5011) != scheduled_nodes.end()) << endl;
                            } else {
                                cout << "Something is very wrong." << endl;
                            }
                            for (auto i = 0; i < num_proc; ++i) {
                                if (end_time_on_proc[i] < min2) {
                                    cout << "Seems like there is a bug." << endl;
                                    exit(1);
                                }
                            }
                        }
                        
                        scheduled_nodes.erase(run_on_proc[min_index2]);
                        sbu_activate_children(run_on_proc[min_index2]);
                        run_on_proc[min_index2] = num_vertices(G) - 1;
                        for (auto j = 0; j < num_proc; ++j) {
                            if (end_time_on_proc[j] == 0) {
                                end_time_on_proc[j] = min2;
                            }
                        }
                        for (auto j = 0; j < num_proc; ++j) {
                            if (end_time_on_proc[j] < min2) {
                                end_time_on_proc[j] = min2;
                                if (run_on_proc[j] != num_vertices(G) - 1) {
                                    cout << "Probably scheduling error (node should had finished)." << endl;
                                    exit(1);
                                }
                            }
                        }
#ifdef DBG_PRINT
                        cout << "===== Current running nodes =====" << endl;
                        for (auto j = 0; j < num_proc; ++j) {
                            cout << "Node ID:\t" << run_on_proc[j] << ".\tEnd time: " << end_time_on_proc[j] << endl;
                        }
#endif
                    } else {
                        break;
                    }
                }
            }
#ifdef DBG_PRINT
            cout << "===== Current running nodes =====" << endl;
            for (auto j = 0; j < num_proc; ++j) {
                cout << "Node ID:\t" << run_on_proc[j] << ".\tEnd time: " << end_time_on_proc[j] << endl;
            }

            cout << "Scheduling node " << *node << " to run..." << endl;
            cout << "Its processing time: " << G[*node].proc_time << endl;
#endif
            // We want to schedule nodes on the processor, which finishes first
            double min = end_time_on_proc[0];
            auto index = 0;
            for (auto i = 1; i < num_proc; ++i) {
                if (end_time_on_proc[i] < min) {
                    min = end_time_on_proc[i];
                    index = i;
                }
            }
#ifdef DBG_PRINT
            cout << "Node " << *node << " scheduled to the processor: " << index << endl;
#endif
            scheduled_nodes.erase(run_on_proc[index]);

//            TODO:
//            - the bug is that it actually go to the next level as all nodes on the
//            current level are scheduled. But they may not complete yet when
//            the new level is being handled
            sbu_activate_children(run_on_proc[index]);
            run_on_proc[index] = *node;
            total_jobs++;
            scheduled_nodes.insert(*node);
            end_time_on_proc[index] += G[*node].proc_time;
            active_nodes_on_level[G[*node].level].erase(*node);
// TODO: check that no jobs activated before jobs has finished
// TODO: check that the last jobs finished and there is nothing more to activate
// TODO: do not skip first level if it empty as may skip doing useful work by allocating all jobs beforehand

#ifdef DBG_PRINT
            cout << "===== Current running nodes =====" << endl;
            for (auto j = 0; j < num_proc; ++j) {
                cout << "Node ID:\t" << run_on_proc[j] << ".\tEnd time: " << end_time_on_proc[j] << endl;
            }
#endif
        }
#ifdef DBG_PRINT
        cout << "===== Finished level " << level << " =====" << endl;
        cout << "===== Current running nodes after finishing level =====" << endl;
        for (auto j = 0; j < num_proc; ++j) {
            cout << "Node ID:\t" << run_on_proc[j] << ".\tEnd time: " << end_time_on_proc[j] << endl;
        }
#endif
    }

    double max = 0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] > max) {
            max = end_time_on_proc[i];
        }
    }

    cout << "Total run time: " << max << endl;
    cout << "Total jobs scheduled: " << total_jobs << endl;

    return ;
}

void run_sbu_blocking_scheduler(uint64_t num_levels) {
#ifdef DEBUG_PRINT
    cout << "Starting SBU blocking scheduler." << endl;
#endif
    auto total_jobs = 0;
    auto active_jobs = 0;

    for (auto i = 1; i < num_levels; ++i) {
        active_jobs += active_nodes_on_level[i].size();
    }
    cout << "Total work span: " << active_jobs << endl;
    sleep(5);

    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc);
    vector<double> end_time_on_proc(num_proc);
    for (auto i = 0; i < num_proc; ++i) {
        run_on_proc[i] = sink_node;
        end_time_on_proc[i] = 0;
    }

    uint64_t init_num_nodes = (active_nodes_on_level[1].size() < num_proc) ? active_nodes_on_level[1].size() : num_proc;
    if (init_num_nodes == 0) {
        cout << "WARNING: level 1 is empty." << endl;
        sleep(3);
    }
    for (auto i = 0; i < init_num_nodes; ++i) {
        auto node = active_nodes_on_level[1].begin();
        run_on_proc[i] = *node;
        end_time_on_proc[i] = G[*node].proc_time;
        active_nodes_on_level[1].erase(*node);
        total_jobs++;
    }

    for (auto level = 1; level < num_levels; ++level) {
        cout << "Processing level " << level << " with " << active_nodes_on_level[level].size() << " active nodes." << endl;
        if (active_nodes_on_level[level].size() == 0) {
            if (level != 1) {
                cout << "Level is empty. Skipping..." << endl;
                continue;
            }
        }
        for (auto node : active_nodes_on_level[level]) {
            if (G[node].is_predicate_node) {
                cout << "ERROR: predicate node." << endl;
                exit(1);
            }
            if (!G[node].is_active) {
                cout << "ERROR: inactive node." << endl;
                exit(1);
            }
        }

        while (!active_nodes_on_level[level].empty()) {
            auto node = active_nodes_on_level[level].begin();

            // We want to schedule nodes on the processor, which finishes first
            double min = end_time_on_proc[0];
            auto index = 0;
            for (auto i = 1; i < num_proc; ++i) {
                if (end_time_on_proc[i] < min) {
                    min = end_time_on_proc[i];
                    index = i;
                }
            }

//            TODO:
//            - the bug is that it actually go to the next level as all nodes on the
//            current level are scheduled. But they may not complete yet when
//            the new level is being handled
            sbu_activate_children(run_on_proc[index]);
            run_on_proc[index] = *node;
            end_time_on_proc[index] += G[*node].proc_time;
            total_jobs++;
            active_nodes_on_level[G[*node].level].erase(*node);
        }

        // As this is blocking algorithm before proceeding to the next level
        // we first should finish all the jobs on the current level. At this point
        // the set of active nodes on the level is empty but there may be still
        // some running nodes
        double max = end_time_on_proc[0];
        auto index = 0;
        for (auto i = 1; i < num_proc; ++i) {
            if (end_time_on_proc[i] > max) {
                max = end_time_on_proc[i];
                index = i;
            }
        }
        for (auto i = 0; i < num_proc; ++i) {
            end_time_on_proc[i] = max;
            sbu_activate_children(run_on_proc[i]);
            run_on_proc[i] = sink_node;
        }
    }

    double max = 0;
    for (auto i = 0; i < num_proc; ++i) {
        if (end_time_on_proc[i] > max) {
            max = end_time_on_proc[i];
        }
    }

    cout << "Total run time: " << max << endl;
    cout << "Total jobs scheduled: " << total_jobs << endl;

    return ;
}

uint64_t run_sbu_blocking_scheduler_measure(uint64_t num_levels) {
    auto total_jobs = 0;
    auto active_jobs = 0;

//    for (auto i = 1; i < num_levels; ++i) {
//        active_jobs += active_nodes_on_level[i].size();
//    }

    vector<DirectedGraph::vertex_descriptor> run_on_proc(num_proc);
    vector<double> end_time_on_proc(num_proc);
    /*
    for (auto i = 0; i < num_proc; ++i) {
        run_on_proc[i] = sink_node;
        end_time_on_proc[i] = 0;
    }
    */

    uint64_t init_num_nodes = (active_nodes_on_level[1].size() < num_proc) ? active_nodes_on_level[1].size() : num_proc;
    /*
    for (auto i = 0; i < init_num_nodes; ++i) {
        auto node = active_nodes_on_level[1].begin();
        run_on_proc[i] = *node;
        end_time_on_proc[i] = G[*node].proc_time;
        active_nodes_on_level[1].erase(*node);
        total_jobs++;
    }
    */

    /*
    auto total_active_jobs_measure = 0;
    for (auto level = 1; level < num_levels; ++level) {
        if (active_nodes_on_level[level].size() == 0) {
            if (level != 1) {
                continue;
            }
        }
        for (auto node : active_nodes_on_level[level]) {
            if (G[node].is_predicate_node) {
                cout << "ERROR: predicate node." << endl;
                exit(1);
            }
            if (!G[node].is_active) {
                cout << "ERROR: inactive node." << endl;
                exit(1);
            }
            total_active_jobs_measure++;
        }
    }
    return total_active_jobs_measure;
    */
    for (auto level = 1; level < num_levels; ++level) {
        if (active_nodes_on_level[level].size() == 0) {
            if (level != 1) {
                continue;
            }
        }
        for (auto node : active_nodes_on_level[level]) {
            if (G[node].is_predicate_node) {
                cout << "ERROR: predicate node." << endl;
                exit(1);
            }
            if (!G[node].is_active) {
                cout << "ERROR: inactive node." << endl;
                exit(1);
            }
        }

        while (!active_nodes_on_level[level].empty()) {
            auto node = active_nodes_on_level[level].begin();

            // We want to schedule nodes on the processor, which finishes first
            double min = end_time_on_proc[0];
            auto index = 0;
            /*
            for (auto i = 1; i < num_proc; ++i) {
                if (end_time_on_proc[i] < min) {
                    min = end_time_on_proc[i];
                    index = i;
                }
            }
            */

//            TODO:
//            - the bug is that it actually go to the next level as all nodes on the
//            current level are scheduled. But they may not complete yet when
//            the new level is being handled
            sbu_activate_children(run_on_proc[index]);
            run_on_proc[index] = *node;
            end_time_on_proc[index] += G[*node].proc_time;
            total_jobs++;
            active_nodes_on_level[G[*node].level].erase(*node);
        }

        // As this is blocking algorithm before proceeding to the next level
        // we first should finish all the jobs on the current level. At this point
        // the set of active nodes on the level is empty but there may be still
        // some running nodes
        double max = end_time_on_proc[0];
        auto index = 0;
        /*
        for (auto i = 1; i < num_proc; ++i) {
            if (end_time_on_proc[i] > max) {
                max = end_time_on_proc[i];
                index = i;
            }
        }
        */
        for (auto i = 0; i < num_proc; ++i) {
            end_time_on_proc[i] = max;
            sbu_activate_children(run_on_proc[i]);
            run_on_proc[i] = sink_node;
        }
    }

    double max = 0;
//    for (auto i = 0; i < num_proc; ++i) {
//        if (end_time_on_proc[i] > max) {
//            max = end_time_on_proc[i];
//        }
//    }

//    cout << "Total run time: " << max << endl;
//    cout << "Total jobs scheduled: " << total_jobs << endl;

    return 0;
}

void sbu_populate_active_list() {
    auto count = 0;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    for (; vit != vend; ++vit) {
        if (G[*vit].is_active) {
            if (*vit == 61135) {
                cout << "found our node" << endl;
                DirectedGraph::in_edge_iterator eit, eend;
                tie(eit, eend) = in_edges(*vit, G);
                for (; eit != eend; ++eit) {
                    cout << "source level: " << G[source(*eit, G)].level << endl;
                    cout << "active? " << G[source(*eit, G)].is_active << " ID: " << source(*eit, G) << endl;
                    cout << "predicate? " << G[source(*eit, G)].is_predicate_node << endl;
                }
                for (auto node : active_nodes_on_level[10]) {
                    if (node == 11483) {
                        cout << "NODE ID: " << node << endl;
                    }
                }
                for (auto node : levels[10]) {
                    if (*node == 11483) {
                        cout << "LEVEL ID: " << *node << endl;
                    }
                }
                sleep(5);
            }
            count++;
            active_nodes_on_level[G[*vit].level].insert(*vit);
        }
    }
#ifdef DEBUG_PRINT
    cout << "Added " << count << " jobs to active list." << endl;
#endif
}

void check_trace() {
    cout << "\033[1;32m===== Analyzing trace file =====\033[0m" << endl;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    auto num_active_jobs = 0;
    auto num_active_pred = 0;
    auto num_active_unit = 0;
    auto num_levels = 0;

    // Check that there are no active nodes on the level zero
    for (auto node : active_nodes_on_level[0]) {
        if (G[node].is_active) {
            cout << "ERROR: found active node on the level zero." << endl;
            exit(1);
        }
    }

    for (; vit != vend; ++vit) {
        if (G[*vit].is_active) {
            num_active_jobs++;
            if (G[*vit].is_predicate_node) {
//                cout << "None of the predicate nodes should be active in the trace file." << endl;
//                exit(1);
                num_active_pred++;
            } else {
                num_active_unit++;
            }
        }
    }
    cout << "Number of active jobs:\t\t\t" << num_active_jobs << "/" << num_active_unit << endl;
    cout << "Number of active predicate nodes:\t" << num_active_pred << "/" << num_active_unit << endl;
    cout << "Number of active unit nodes:\t\t" << num_active_unit << endl;

    num_levels = find_depth();
    auto active_through_levels = 0;
    for (auto i = 0; i < num_levels; ++i) {
        for (auto el : levels[i]) {
            if (G[*el].is_active) {
                active_through_levels++;
            }
        }
    }
    cout << "Number of active jobs (iterating through levels): " << active_through_levels << endl;

    for (auto el : levels[0]) {
        if (G[*el].is_scheduled) {
            cout << "ERROR: There is a scheduled node on level zero." << endl;
            exit(1);
        }
    }

    for (auto i = 0; i < num_levels; ++i) {
        for (auto el : levels[i]) {
            if (G[*el].is_scheduled) {
                DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
                tie(edgeIt, edgeEnd) = out_edges(*el, G);
                for (; edgeIt != edgeEnd; ++edgeIt) {
                    if (!G[target(*edgeIt, G)].is_active) {
                        cout << "There is a node of a scheduled node, which is not active." << endl;
                        exit(1);
                    }
                }
            }
        }
    }
    cout << "All of the children nodes of scheduled nodes are active." << endl;

    cout << "\033[1;32m================================\033[0m" << endl;
}

void activate_initial_nodes() {
    auto num_levels = find_depth();
    bfs_active_nodes = 0;
    cout << "Number of nodes on level 1: " << levels[1].size() << endl;
    uint64_t num_init_scheduled = 0;
    for (auto el : levels[1]) {
        if (G[*el].is_scheduled) {
            DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
            tie(edgeIt, edgeEnd) = out_edges(*el, G);
            for (; edgeIt != edgeEnd; ++edgeIt) {
                if (!G[target(*edgeIt, G)].is_active) {
                    cout << "Trying to schedule not active node." << endl;
                    exit(1);
                }
                G[target(*edgeIt, G)].is_scheduled = true;
                num_init_scheduled++;
                scheduled_nodes.insert(target(*edgeIt, G));
            }
            counting_bfs vis;
            breadth_first_search(G, *el, visitor(vis));
        }
        if (G[*el].is_active && !G[*el].is_predicate_node) {
            cout << "ERROR: Active unit node on level 1." << endl;
            exit(1);
        }
    }
#ifdef DEBUG_PRINT
    cout << "Number of initial nodes scheduled:" << num_init_scheduled << endl;
    cout << "Active nodes reachable from scheduled source nodes (with repetitions): " << bfs_active_nodes <<  endl;
    cout << "Size of the reachable active set: " << bfs_active_set.size() << endl;
#endif

}

double critical_path_routine(DirectedGraph::vertex_descriptor v) {
    DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
    tie(edgeIt, edgeEnd) = out_edges(v, G);
    double max = 0;
    for (; edgeIt != edgeEnd; ++edgeIt) {
        if (!G[target(*edgeIt, G)].is_active) {
            continue;
        }
        double path_cost = critical_path_routine(target(*edgeIt, G));
        if (path_cost > max) {
            max = path_cost;
        }
    }
    if (max != 0) {
        return G[v].proc_time + max;
    } else {
        return G[v].proc_time;
    }
}

uint64_t num_active_nodes_on_level(uint64_t level) {
    uint64_t count = 0;
    for (auto el : levels[level]) {
        if (G[*el].is_active) {
            count++;
        }
    }

    return count;
}

uint64_t num_scheduled_nodes_on_level(uint64_t level) {
    uint64_t count = 0;
    for (auto el : levels[level]) {
        if (G[*el].is_scheduled) {
            count++;
        }
    }

    return count;
}

void find_critical_path() {
    double max = 0;
    uint64_t count = 0;
    for (auto el: levels[1]) {
        if (!G[*el].is_scheduled) {
            continue;
        }
        DirectedGraph::out_edge_iterator edgeIt, edgeEnd;
        tie(edgeIt, edgeEnd) = out_edges(*el, G);
        for (; edgeIt != edgeEnd; ++edgeIt) {
            double path_cost = critical_path_routine(target(*edgeIt, G));
            cout << "Path cost: " << path_cost << endl;
            if (path_cost > max) {
                max = path_cost;
            }
        }
        cout << "Critical path: Checked " << ++count << "/" << num_scheduled_nodes_on_level(1) << ".\r" << flush;
    }
    cout << endl << "Optimal running time with infinite parallelism: " << max << " seconds." << endl;
}

void check_graph(vector<DirectedGraph::vertex_iterator> source_nodes) {
    cout << "Performing integrity checks..." << endl;
    // TODO:
    // - check number of predicate nodes in the graph
    // - check number of ingress/egress edges for predicate and unit nodes
    // - check children of predicate and unit nodes
    // - check that all active nodes are reachable (code snippet is commented below)
    // - check that there are no active source unit nodes
    // - here should be all checks performed
    // - check that there are no active predicate nodes
    DirectedGraph::vertex_iterator vertexIt, vertexEnd;
    tie(vertexIt, vertexEnd) = vertices(G);
    for (; vertexIt != vertexEnd; ++vertexIt) {
        if (G[*vertexIt].is_predicate_node) {
            if (in_degree(*vertexIt, G) > 1) {
                cout << "ERROR: Predicate node has in-degree of more than one." << endl;
            }
        }
    }
    cout << "check source nodes" << endl;
    for (auto node : source_nodes) {
        if (!G[*node].is_predicate_node) {
            DirectedGraph::out_edge_iterator eit, eend;
            tie(eit, eend) = out_edges(*node, G);
            for ( ; eit != eend; ++eit) {
                if (!G[target(*eit, G)].is_predicate_node) {
                    cout << "httt" << endl;
                }
            }
        }
    }
    /*
    for (auto i = 2; i < num_levels; ++i) {
        for (auto el : levels[i]) {
            if (G[*el].is_active) {
                set<DirectedGraph::vertex_descriptor>::iterator it;
                it = bfs_active_set.find(*el);
                if (it == bfs_active_set.end()) {
                    cout << "ERROR: There is a non-reachable active node." << endl;
                    exit(1);
                }
            }
        }
    }
    */
    cout << "edges: " << num_edges(G) << endl;
    cout << "Vertices: " << num_vertices(G) << endl;
}

void check_reachability(vector<DirectedGraph::vertex_iterator> source_nodes) {
    size_t num_levels = find_depth();
    for (auto i = 1; i < num_levels; ++i) {
        for (auto node : levels[i]) {
            if (G[*node].is_scheduled) {
                if (G[*node].prunned) {
                    cout << "ERROR: pruned scheduled node." << endl;
                    exit(1);
                }
                sbu_activate_children(*node);
            }
        }
    }

    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(G);
    for ( ; vit != vend; ++vit) {
        if (!G[*vit].prunned) {
            if (G[*vit].level == 0) {
                continue;
            }
            if (G[*vit].is_active && !G[*vit].is_scheduled) {
                cout << "Active but not scheduled node with ID: " << *vit << endl;
                cout << "Its level: " << G[*vit].level << endl;
                exit(1);
            }
        }
    }

    cout << "All active nodes are reachable." << endl;
    exit(1);
}

void extract_subgraph(DirectedGraph graph, vector<DirectedGraph::vertex_iterator> source_nodes) {
    subgraph_bfs vis;
    vector<DirectedGraph::vertex_iterator> tmp_subgraph_nodes;
    for (auto node : source_nodes) {
        if (*node == 36345 || *node == 34743 || *node == 54196 || *node == 27376 || *node == 3230) {
            tmp_subgraph_nodes.push_back(node);
            breadth_first_search(graph, *node, visitor(vis));
        }
    }
    cout << "Checking correctness: " << tmp_subgraph_nodes.size() << "/" << source_nodes.size() << endl;
    DirectedGraph::vertex_iterator vit, vend;
    tie(vit, vend) = vertices(graph);
    for ( ; vit != vend; ++vit) {
        auto tmp = subgraph_nodes.find(*vit);
        if (tmp != subgraph_nodes.end()) {
            tmp_subgraph_nodes.push_back(vit);
        }
    }
    vertices_count = tmp_subgraph_nodes.size();
//    subgraph<DirectedGraph>& mygraph = graph.create_subgraph(tmp_subgraph_nodes.begin(), tmp_subgraph_nodes.end());
//    cout << "Subgraph has " << num_vertices(mygraph) << " vertices and " << num_edges(mygraph) << " edges." << endl;
    cout << "Subgraph has " << vertices_count << " vertices and " << edges_count << " edges." << endl;
    edges_count = subgraph_edges.size();
    vertices_count = subgraph_nodes.size();
    cout << "Subgraph has " << vertices_count << " vertices and " << edges_count << " edges." << endl;
    auto pnode = 0;
    auto unode = 0;
    for (auto node : subgraph_nodes) {
        if (graph[node].is_predicate_node) {
            pnode++;
        } else {
            unode++;
        }
    }
    cout << "predicate nodes: " << pnode << endl << "unit nodes: " << unode << endl;
//    exit(1);
}

void finish_lookahead_earliest_job(vector<double> &end_time_on_proc, vector<DirectedGraph::vertex_descriptor> &run_on_proc, bool &CPUfree, bool &allCPUfree) {
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

    sbu_activate_children(run_on_proc[cpu_num]);
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

    run_on_proc[cpu_num] = sink_node;
    CPUfree = true;
    num_executed_jobs++;

    for (auto i = 0; i < num_proc; ++i) {
        if (run_on_proc[i] != sink_node)
            return ;
    }
    allCPUfree = true;
}

DirectedGraph::vertex_descriptor get_lookahead_scheduled_node() {
    sched_start = chrono::high_resolution_clock::now();
    DirectedGraph::vertex_descriptor node = -1;

    // First run nodes on the current level
    if (!run_queue.empty()) {
        node = *(run_queue.begin());
//        run_queue.erase(node);

        sched_end = chrono::high_resolution_clock::now();
        scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

        return node;
    }

    // If there are no nodes on the current level then run lookahead
    uint64_t tmp_level2 = current_level + 1;
    for (auto tmp_level = tmp_level2; tmp_level < tmp_level2 + lookahead_levels; tmp_level++) {
        if (tmp_level >= find_depth())
            return -1;
        for (auto l_node = active_nodes_on_level[tmp_level].begin(); l_node != active_nodes_on_level[tmp_level].end(); ++l_node) {
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
                return *l_node;
            }
        }
    }

    sched_end = chrono::high_resolution_clock::now();
    scheduling_overhead += chrono::duration_cast<chrono::microseconds>(sched_end - sched_start);

    return node;
}

void run_sbu_lookahead_scheduler(size_t num_levels) {
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
        auto node = get_lookahead_scheduled_node();
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
            finish_lookahead_earliest_job(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
        }
        auto node = get_lookahead_scheduled_node();
        if (scheduled_nodes.find(node) == scheduled_nodes.end() && node != -1) {
            cout << "ERROR: Trying to schedule node not from the scheduled set (2)." << endl;
            exit(1);
        }
        while (node == -1) {
            cout << "Trying to schedule by finishing work..." << endl;
            finish_lookahead_earliest_job(end_time_on_proc, run_on_proc, CPUfree, allCPUfree);
            cout << "Finished the job." << endl;
            node = get_lookahead_scheduled_node();
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

void sbu_scheduler(bool trace_enabled) {
    cout << "Starting SBU scheduler..." << endl;
    cout << "vertices: " << num_vertices(G) << endl;
    cout << "edges: " << num_edges(G) << endl;
    vector<DirectedGraph::vertex_iterator> source_nodes;
    sbu_init_graph(source_nodes);
    check_graph(source_nodes);
//        extract_subgraph(G, source_nodes);
    rebuild_graph();
    find_levels_ds(G_tmp);

    size_t num_levels = find_depth();
    levels.resize(num_levels + 1);
    populate_levels();
    running_nodes_on_level.resize(num_levels + 1);
    active_nodes_on_level.resize(num_levels + 1);

    if (trace_enabled) {
        add_start_nodes();


        assign_jobs_from_trace();
        cout << "Number of look ahead levels: " << lookahead_levels << endl;
        cout << "vertices: " << num_vertices(G) << endl;
        cout << "edges: " << num_edges(G) << endl;
        cout << "vertices: " << num_vertices(G_tmp) << endl;
        cout << "edges: " << num_edges(G_tmp) << endl;
        cout << "NUMBER OF SOURCE NODES: " << source_nodes.size() << endl;
        cout << "NUMBER OF LEVELS NODES: " << levels[0].size() << endl;
        for (auto node : levels[0]) {
            if (in_degree(*node, G) != 0) {
                cout << "IN DEGREE NOT ZERO." << endl;
            }
        }
        check_trace(); // TODO: check this function
    } else {
        assign_jobs();
    }

    sbu_populate_active_list();

    add_sink_node();
    num_levels = find_depth();
    /*
    sleep(5);
    if (draw_graph(0, levels, num_levels, "graph.gv") < 0) {
        cout << "Error during producing file for graph plotting." << endl;
        exit(0);
    }
    cout << "Done drawing." << endl;
    exit(0);
    */

    cout << "Starting scheduler..." << endl;
    for (auto node : source_nodes) {
        sbu_activate_children(*node);
    }

    //check_reachability(source_nodes);

    /*
    bool found_single_node = false;
    for (auto node : source_nodes) {
        node_in_search = 61124;
        node_in_search_found = false;
        searching_bfs vis;
        breadth_first_search(G, *node, visitor(vis));
        if (node_in_search_found) {
            cout << "Node is reachable." << endl;
            found_single_node = true;
            break;
        }
    }
    if (!found_single_node) {
        cout << "ERROR: The node in question is not reachable from source nodes." << endl;
        exit(1);
    }
    */

//    cout << "Searching for critical path..." << endl;
//    find_critical_path();
//    exit(1);
//    fully_parallel_experiment(num_levels + 1);
//    run_sbu_scheduler(num_levels + 1);
    run_sbu_lookahead_scheduler(num_levels + 1);
    
//        run_sbu_blocking_scheduler(num_levels + 1);
//    auto start = chrono::steady_clock::now();
//    auto res = run_sbu_blocking_scheduler_measure(num_levels + 1);
//    auto end = chrono::steady_clock::now();
//    cout << "Number of touched active jobs: " << res << endl;

    /*
    cout << "Elapsed time in microseconds : " 
        << chrono::duration_cast<chrono::microseconds>(end - start).count()
        << " µs" << endl;

    cout << "Elapsed time in milliseconds : " 
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << endl;

    cout << "Elapsed time in seconds : " 
        << chrono::duration_cast<chrono::seconds>(end - start).count()
        << " sec" << endl;
    */
    
    //    elapsed_time(start, end);
//    find_critical_path();
#ifdef DEBUG_PRINT
    cout << "Approximate time in units to run is: " << calculate_time(levels, num_levels + 1) << endl;
#endif
}
