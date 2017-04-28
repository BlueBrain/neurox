COMMENT
/**
 * @file MemUsage.mod
 * @brief 
 * @author king
 * @date 2011-02-04
 * @remark Copyright © BBP/EPFL 2005-2011; All rights reserved. Do not distribute without further notice.
 */
ENDCOMMENT

VERBATIM
#include <mpi.h>
#include <malloc.h> //for mallinfo
ENDVERBATIM

NEURON {
    ARTIFICIAL_CELL MemUsage
}


PARAMETER {
    minUsageMB = 0
    maxUsageMB = 0
    avgUsageMB = 0
    stdevUsageMB = 0
    rank = 0
    size = 0
}

INITIAL {
}

CONSTRUCTOR  {
VERBATIM
    int i_rank, i_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &i_size);    
    rank = i_rank;
    size = i_size;
ENDVERBATIM
}

PROCEDURE print_mem_usage() {
VERBATIM 
/**
 * Gather memory usage statistics for all nodes in the network, printing to the console
 */
    struct mallinfo memInfo = mallinfo();
    double usageMB = (double) (memInfo.arena + memInfo.hblkhd)/ (double) (1024*1024);

    MPI_Reduce( &usageMB, &minUsageMB, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &usageMB, &maxUsageMB, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &usageMB, &avgUsageMB, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    avgUsageMB /= size;
    
    MPI_Bcast( &avgUsageMB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
    double diffSquared = (usageMB-avgUsageMB)*(usageMB-avgUsageMB);
    MPI_Reduce( &diffSquared, &stdevUsageMB, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    stdevUsageMB = sqrt( stdevUsageMB/size);
    
    if( rank == 0 ) {
        printf( " memusage according to mallinfo:\n\tMax = %lfMB\n\tMin = %lfMB\n\tMean(Stdev) = %lfMB(%lfMB)\n\n", maxUsageMB, minUsageMB, avgUsageMB, stdevUsageMB );
    }
ENDVERBATIM
}

PROCEDURE print_node_mem_usage() {
VERBATIM 
/**
 * Print memory usage statistics for the local node
 */
    struct mallinfo memInfo = mallinfo();
    double usageMB = (double) (memInfo.arena + memInfo.hblkhd)/ (double) (1024*1024);
    printf( " memusage node %.0lf according to mallinfo:\n\t %lfMB\n\n", rank, usageMB );
ENDVERBATIM
}
