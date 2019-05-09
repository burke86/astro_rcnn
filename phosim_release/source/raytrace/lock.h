///
/// @package phosim
/// @file lock.h
/// @brief header for lock
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

struct Lock {
    /* pthread_mutex_t lock1; */
    pthread_mutex_t lock2;
    /* pthread_mutex_t lock3; */
    pthread_mutex_t lock4;
    pthread_mutex_t lock5;
    pthread_mutex_t lock6;
    pthread_mutex_t lock7;
    pthread_mutex_t lock8;
    pthread_cond_t cond;
    /* pthread_cond_t cond2; */
};
