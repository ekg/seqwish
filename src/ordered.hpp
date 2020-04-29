#pragma once

#include <thread>
#include <atomic>
#include <cstdint>
#include <functional>
#include <map>
#include <utility>
#include "atomic_queue.h"

namespace seqwish {

template<typename Order, typename Value>
void parallel_ordered(
    const Order& begin,
    const Order& end,
    const uint64_t& num_threads,
    const std::function<std::pair<Order, Value*>(const Order&)>& process,
    const std::function<void(const std::pair<Order, Value*>&)>& emit) {
    
    // producer/consumer queues
    auto todo_q_ptr = new atomic_queue::AtomicQueue2<Order, 2 << 16>;
    auto& todo_q = *todo_q_ptr;
    auto done_q_ptr = new atomic_queue::AtomicQueue2<std::pair<Order, Value*>, 2 << 16>;
    auto& done_q = *done_q_ptr;
    std::atomic<bool> work_todo;
    std::map<Order, Value*> node_records;

    auto worker_lambda =
        [&](void) {
            Order id = 0;
            while (work_todo.load()) {
                if (todo_q.try_pop(id)) {
                    done_q.push(process(id));
                } else {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                }
            }
        };

    std::vector<std::thread> workers; workers.reserve(num_threads);
    work_todo.store(true);
    for (uint64_t t = 0; t < num_threads; ++t) {
        workers.emplace_back(worker_lambda);
    }
    Order todo_id = begin;
    Order done_id = 0;
    while (done_id < end) {
        // put ids in todo
        while (todo_id <= end && todo_q.try_push(todo_id)) {
            ++todo_id;
        }
        // read from done queue
        std::pair<Order, Value*> item;
        while (done_q.try_pop(item)) {
            node_records[item.first] = item.second;
        }
        if (node_records.size()) {
            auto b = node_records.begin();
            if (b->first == done_id+1) {
                emit(*b);
                ++done_id;
                delete b->second;
                node_records.erase(b);
            }
        }
    }
    work_todo.store(false);
    for (uint64_t t = 0; t < num_threads; ++t) {
        workers[t].join();
    }
}

}
