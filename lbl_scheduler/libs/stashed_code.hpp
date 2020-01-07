void exec_job(int level) {
    this_thread::sleep_for(chrono::milliseconds(10));
    lock_guard<mutex> lock(jobs_mutex);
    auto node = running_nodes_on_level[level].begin();
    activate_children(*node);
    running_nodes_on_level[level].erase(*node);
#ifdef DEBUG_PRINT
    cout << "Job has been finished on the thread " << std::this_thread::get_id() << endl;
    cout << "Number of running nodes on the current level is: " << running_nodes_on_level[level].size() << endl;
    cout << "Number of active jobs on the current level is: " << active_nodes_on_level[level].size() << endl;
#endif
}

void exec_job_single_thread(int level) {
    this_thread::sleep_for(chrono::milliseconds(10));
    // lock_guard<mutex> lock(jobs_mutex);
    auto node = running_nodes_on_level[level].begin();
    activate_children(*node);
    running_nodes_on_level[level].erase(*node);
#ifdef DEBUG_PRINT
    cout << "Job has been finished on the thread " << std::this_thread::get_id() << endl;
    cout << "Number of running nodes on the current level is: " << running_nodes_on_level[level].size() << endl;
    cout << "Number of active jobs on the current level is: " << active_nodes_on_level[level].size() << endl;
#endif
}

void run(set<uint64_t> &active_jobs, int level) {
    // lock_guard<mutex> lock(jobs_mutex);
    auto node = active_jobs.begin();
    running_nodes_on_level[level].insert(*node);
    active_jobs.erase(*node);
    exec_job_single_thread(level);
    // std::thread thr = thread(&exec_job, level);
    // thr.detach();
}

void scheduler(uint64_t num_levels) {
#ifdef DEBUG_PRINT
    cout << "Starting scheduler()." << endl;
#endif

    running_nodes_on_level.resize(num_levels);

    for (auto i = 1; i < num_levels; ++i) {
#ifdef DEBUG_PRINT
        cout << "Starting to process level " << i << endl;
        cout << "Number of active nodes on the current level is: " << active_nodes_on_level[i].size() << endl;
        cout << "Number of running nodes on the current level is: " << running_nodes_on_level[i].size() << endl;
#endif
        while(!active_nodes_on_level[i].empty() || !running_nodes_on_level[i].empty()) {
            // cout << "Are active nodes empty: " << active_nodes_on_level[i].empty() << endl;
            // cout << "Are running nodes empty: " << running_nodes_on_level[i].empty() << endl;
            if (!active_nodes_on_level[i].empty()) {
#ifdef DEBUG_PRINT
                cout << "Active nodes queue is not empty. Trying to schedule the job (number of active jobs on the current level is "
                     << active_nodes_on_level[i].size() << ")." << endl;
                cout << "Number of running nodes on the current level is: " << running_nodes_on_level[i].size() << endl;
#endif
                run(active_nodes_on_level[i], i);
            }
        }
#ifdef DEBUG_PRINT
        cout << "Level " << i << " is complete." << endl;
#endif
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
            output_file << *vertexIt << endl;
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

void run_job(int i, int level, vector<bool> &thread_running) {
    this_thread::sleep_for(chrono::milliseconds(10));
    lock_guard<mutex> lock(jobs_mutex);
    // running_nodes_on_level[level].pop_back();
    auto node = running_nodes_on_level[level].begin();
    // activate_children(*node);
    running_nodes_on_level[level].erase(*node);
    thread_running[i] = false;
#ifdef DEBUG_PRINT
    cout << "Job has been finished on the thread " << i << endl;
#endif
}

void schedule_job(set<uint64_t> &active_jobs, const int &num_threads, thread threads[], const int &level, vector<bool> &thread_running) {
#ifdef DEBUG_PRINT
    cout << "Searching for the available thread to schedule the job..." << endl;
#endif
    for (auto i = 0; i < num_threads; ++i) {
        if (!thread_running[i]) {
#ifdef DEBUG_PRINT
            cout << "Thread " << i << " is not running. Scheduling job on it..." << endl;
#endif
            thread_running[i] = true;
            {
                lock_guard<mutex> lock(jobs_mutex);
                // running_nodes_on_level[level].push_back(active_jobs.back());
                auto node = active_jobs.begin();
                running_nodes_on_level[level].insert(*node);
                active_jobs.erase(*node);
            }
//            threads[i] = thread(&run_job, i, level, std::ref(thread_running));
//            threads[i].detach();
#ifdef DEBUG_PRINT
            cout << "Job was scheduled on the thread " << i << endl;
#endif
            return ;
        }
    }
#ifdef DEBUG_PRINT
    cout << "No available threads to run jobs." << endl;
#endif
}

void run_old_scheduler(uint64_t num_levels) {
#ifdef DEBUG_PRINT
    cout << "Starting SBU scheduler." << endl;
#endif
    int num_threads = std::thread::hardware_concurrency();
    if (num_threads < 2) {
        cout << "WARNING: Target machine have less than 2 cores." << endl;
    }
    vector<bool> thread_running(num_threads - 1);
    thread *threads = new thread[num_threads - 1];
    running_nodes_on_level.resize(num_levels);
    for (auto i = 0; i < num_threads - 1; ++i) {
        thread_running[i] = false;
    }
#ifdef DEBUG_PRINT
    cout << "Number of threads: " << num_threads << endl;
#endif

    for (auto i = 1; i < num_levels; ++i) {
#ifdef DEBUG_PRINT
        cout << "Starting to process level " << i << endl;
        cout << "Number of active nodes on the current level is: " << active_nodes_on_level[i].size() << endl;
        cout << "Number of running nodes on the current level is: " << running_nodes_on_level[i].size() << endl;
#endif
        while(!active_nodes_on_level[i].empty() || !running_nodes_on_level[i].empty()) {
            if (!active_nodes_on_level[i].empty()) {
#ifdef DEBUG_PRINT
                cout << "Active nodes queue is not empty. Trying to schedule the job (number of active jobs on the current level is "
                     << active_nodes_on_level[i].size() << ")." << endl;
                cout << "Number of running nodes on the current level is: " << running_nodes_on_level[i].size() << endl;
#endif
                int level = i;
                schedule_job(active_nodes_on_level[i], num_threads - 1, threads, level, thread_running);
            }
        }
#ifdef DEBUG_PRINT
        cout << "Level " << i << " is complete." << endl;
#endif
    }

    delete[] threads;
}


