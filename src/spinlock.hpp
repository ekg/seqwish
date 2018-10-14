#ifndef SPINLOCK_HPP_
#define SPINLOCK_HPP_

#include <atomic>

namespace seqwish {

class SpinLock {
public:
    void lock(void) { while(lck.test_and_set(std::memory_order_acquire)) { } }
    void unlock(void) { lck.clear(std::memory_order_release); }
private:
    std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

}

#endif /* SPINLOCK_HPP_ */
