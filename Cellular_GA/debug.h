/**
 * Printing debugging information
 */

#ifndef __DEBUG_H__
#define __DEBUG_H__

#define __MY_DEBUG_VOID_CAST__ (void)

#ifdef NDEBUG_ADDITIONAL
// We don't want any debugging
#  define dbg_vmsg(a)       __MY_DEBUG_VOID_CAST__  0

#  define dbg_msg(a)        __MY_DEBUG_VOID_CAST__  0

#  define dbg_func(a)       __MY_DEBUG_VOID_CAST__  0

#else
    // Perform debugging
#   include <iostream>

    /*********************************************************************
    * prints debugging message to stderr
    */
#ifdef ARC_MPI_ENABLED
#   define dbg_msg(a)                                                                       \
    {                                                                                       \
        int __rank;                                                                         \
        MPI_Comm_rank(MPI_COMM_WORLD, &__rank);                                             \
        std::cerr << "Process " << __rank << ": " <<  __func__ << ": " << a << std::flush;  \
    }
#else
#   define dbg_msg(a)      std::cerr << __func__ << ": " << a << std::flush;
#endif /* ARC_MPI_ENABLED */
#   define dbg_func(a) a

#endif

#endif // __DEBUG_H__

// End of file debug.h
