#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int min(int a, int b) {
    return (a > b) ? b: a;
}

int main(int argc,char* argv[]) {


if (argc != 3) {
        perror("./tema3 N ERROR\n");
    }
    int rank, procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int topologie[3][procs - 3];
    for(int i = 0; i < 3; i++) {
        for (int j = 0; j < procs - 3; j++) {
            topologie[i][j] = 0;
        }
    }
    int *vec0;
    int *vec1;
    int *vec2;
    int nr_proc1 = 0, nr_proc2 = 0, nr_proc3 = 0;
    int N;
    int *calculateVector;
    int start, end;
    int vecStart[procs - 3];
    int vecEnd[procs - 3];
    int vecSCLS[3];
    int vecECLS[3];

    if (rank == 0) {
        FILE *clst0 = fopen("cluster0.txt", "r");
        fscanf(clst0, "%d", &nr_proc1);
        vec0 = malloc(nr_proc1 * sizeof(int));
        for (int i = 0; i < nr_proc1; i++) {
            fscanf(clst0, "%d", &vec0[i]);
            topologie[rank][i] = vec0[i];
        }
        N = atoi(argv[1]);
        calculateVector = malloc(N * sizeof(int));
        for (int i = 0; i < N; i++) {
            calculateVector[i] = i;
        }
        for (int i = 0; i < procs - 3; i++) {
            vecStart[i] = (i * N) / (procs - 3);
            vecEnd[i] = min((((i + 1) * N) / (procs - 3)), N);
        }
        
        // trimite nr de procese la coord 1
        MPI_Send(&nr_proc1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        // trimite nr de procese la coord 2
        MPI_Send(&nr_proc1, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
        // trimite workerii la coord 1
        MPI_Send(vec0, nr_proc1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        // trimite workerii la coord 2
        MPI_Send(vec0, nr_proc1, MPI_INT, 2, 0, MPI_COMM_WORLD);

        // primeste nr de proc de la coord 1, aloca mem si primeste si workerii
        MPI_Recv(&nr_proc2, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(1,%d)\n", rank);
        vec1 = malloc(nr_proc2 * sizeof(int));
        MPI_Recv(vec1, nr_proc2, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(1,%d)\n", rank);
        for (int i = 0; i < nr_proc2; i++) {
            topologie[1][i] = vec1[i];
        }
        // primeste nr de proc de la coord 2, aloca mem si primeste si workerii
        MPI_Recv(&nr_proc3, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(2,%d)\n", rank);
        vec2 = malloc(nr_proc3 * sizeof(int));
        MPI_Recv(vec2, nr_proc3, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(2,%d)\n", rank);
         for (int i = 0; i < nr_proc3; i++) {
            topologie[2][i] = vec2[i];
        }
        // trimit numarul de numere si vectorul la 1
        MPI_Send(&N, 1, MPI_INT, 1, 0 ,MPI_COMM_WORLD);
        MPI_Send(calculateVector, N, MPI_INT, 1, 0 ,MPI_COMM_WORLD);
        
        // trimit numarul de numere si vectorul la 1
        MPI_Send(&N, 1, MPI_INT, 2, 0 ,MPI_COMM_WORLD);
        MPI_Send(calculateVector, N, MPI_INT, 2, 0 ,MPI_COMM_WORLD);

        // formez start si end pentru clustere si le trimit
        vecSCLS[0] = 0;
        vecECLS[0] = vecEnd[nr_proc1 - 1];
        vecSCLS[1] = vecStart[nr_proc1];
        vecECLS[1] = vecEnd[nr_proc1 + nr_proc2 - 1];
        vecSCLS[2] = vecStart[nr_proc1 + nr_proc2];
        vecECLS[2] = vecEnd[nr_proc1 + nr_proc2 + nr_proc3 - 1];
        MPI_Send(vecSCLS, 3, MPI_INT, 1, 0 ,MPI_COMM_WORLD);
        MPI_Send(vecECLS, 3, MPI_INT, 1, 0 ,MPI_COMM_WORLD);
        MPI_Send(vecSCLS, 3, MPI_INT, 2, 0 ,MPI_COMM_WORLD);
        MPI_Send(vecECLS, 3, MPI_INT, 2, 0 ,MPI_COMM_WORLD);


        // pentru toate procesele din coord 0 trimit toti workerii
        for (int i = 0; i < nr_proc1; i++) {
            MPI_Send(&nr_proc1, 1, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            MPI_Send(&nr_proc2, 1, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            MPI_Send(&nr_proc3, 1, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            for (int j = 0; j < 3; j++) {
                MPI_Send(topologie[j], procs - 3, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            }
            MPI_Send(&vecStart[i], 1, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            MPI_Send(&vecEnd[i], 1, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);
            MPI_Send(calculateVector, N, MPI_INT, vec0[i], 0, MPI_COMM_WORLD);

        }
        // primesc de la fiecare coordonator un vector modificat si il imbin 
        int auxV[N];
        int auxV0[N];
        int auxV1[N];
        int auxV2[N];
        for (int i = 0; i < N ; i++) {
            auxV0[i] = calculateVector[i];
        }    

        for (int i = 0; i < nr_proc1; i++) {
            MPI_Recv(auxV, N, MPI_INT, vec0[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("M(%d,%d)\n",vec0[i], rank);    
            for (int j = vecStart[i]; j < vecEnd[i]; j++) {
                auxV0[j] = auxV[j];
            }
        }
        MPI_Recv(auxV1, N, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(1,%d)\n", rank);
        MPI_Recv(auxV2, N, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(2,%d)\n", rank);

        for (int i = 0; i < N; i++) {
            if (calculateVector[i] != auxV0[i]) {
                calculateVector[i] = auxV0[i];
            } else { 
                if (calculateVector[i] != auxV1[i]) {
                    calculateVector[i] = auxV1[i];
                }else {
                if (calculateVector[i] != auxV2[i]) {
                    calculateVector[i] = auxV2[i];
                    }
                }
            } 
        }
        //printez vectorul 
        printf("Rezultat: ");
        for (int i = 0; i < N; i++) {
            printf("%d ", calculateVector[i]);
        }
        printf("\n");

        // printez topologia
        printf("%d -> ", rank);
        printf("0:");
        for (int j = 0; j < nr_proc1 - 1; j++) {
            printf("%d,", topologie[0][j]);
        }
        printf("%d 1:", topologie[0][nr_proc1 - 1]);
        for (int j = 0; j < nr_proc2 - 1; j++) {
            printf("%d,", topologie[1][j]);
        }
        printf("%d 2:", topologie[1][nr_proc2 - 1]);
        for (int j = 0; j < nr_proc3 - 1; j++) {
            printf("%d,", topologie[2][j]);
        }
        printf("%d", topologie[2][nr_proc3 - 1]);
        printf("\n");
    }   
    if (rank == 1) {
        FILE *clst1 = fopen("cluster1.txt", "r");
        fscanf(clst1, "%d", &nr_proc2);
        vec1 = malloc(nr_proc2 * sizeof(int));
        for (int i = 0; i < nr_proc2; i++) {
            fscanf(clst1, "%d", &vec1[i]);
            topologie[rank][i] = vec1[i];
        }
        // trimite nr de procese la coord 0
        MPI_Send(&nr_proc2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        // trimite nr de procese la coord 2
        MPI_Send(&nr_proc2, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);

        // trimite workerii la coord 0
        MPI_Send(vec1, nr_proc2, MPI_INT, 0, 0, MPI_COMM_WORLD);
        // trimite workerii la coord 2
        MPI_Send(vec1, nr_proc2, MPI_INT, 2, 0, MPI_COMM_WORLD);
      
        // primeste nr de proc de la coord 0, aloca mem si primeste si workerii
        MPI_Recv(&nr_proc1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        vec0 = malloc(nr_proc1 * sizeof(int));
        MPI_Recv(vec0, nr_proc1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        for (int i = 0; i < nr_proc1; i++) {
            topologie[0][i] = vec0[i];
        }
        // primeste nr de proc de la coord 2, aloca mem si primeste si workerii
        MPI_Recv(&nr_proc3, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(2,%d)\n", rank);
        vec2 = malloc(nr_proc3 * sizeof(int));
        MPI_Recv(vec2, nr_proc3, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(2,%d)\n", rank);
        for (int i = 0; i < nr_proc3; i++) {
            topologie[2][i] = vec2[i];
        }
        // primesc numarul si vectorul de calculat
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        calculateVector = malloc(N * sizeof(int));
        MPI_Recv(calculateVector, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);

        // primesc vectorii de start si end ai clusterelor
        MPI_Recv(vecSCLS, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        MPI_Recv(vecECLS, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        // imi setez dimensiunea din vector pentru clusterul 1
        int aux = vecECLS[1] - vecSCLS[1];
        for (int i = 0; i < nr_proc2; i++) {
        vecStart[i] = (i *  aux / nr_proc2) + vecSCLS[1];
        vecEnd[i] = min(((i + 1) * aux) / nr_proc2, aux)  + vecSCLS[1];
        }

        // pentru toate procesele din coord 0 trimit toti workerii
        for (int i = 0; i < nr_proc2; i++) {
            MPI_Send(&nr_proc1, 1, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            MPI_Send(&nr_proc2, 1, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            MPI_Send(&nr_proc3, 1, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            for (int j = 0; j < 3; j++) {
                MPI_Send(topologie[j], procs - 3, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            }
            // trimit valorile de start si de end pentru fiecare worker 
            MPI_Send(&vecStart[i], 1, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            MPI_Send(&vecEnd[i], 1, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);
            MPI_Send(calculateVector, N, MPI_INT, vec1[i], 0, MPI_COMM_WORLD);


        }
        // primesc vectorul schimbat de workeri intr o variabila aux si il modific pe cel original
        // doar in intervalul descoperit mai sus. Trimit noul vector 
        int auxV[N];
        for (int i = 0; i < nr_proc2; i++) {
            MPI_Recv(auxV, N, MPI_INT, vec1[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("M(%d,%d)\n",vec1[i], rank);
            for (int j = vecStart[i]; j < vecEnd[i]; j++) {
                calculateVector[j] = auxV[j];
            }
        }
        MPI_Send(calculateVector, N, MPI_INT, 0, 0, MPI_COMM_WORLD);

        //printez topologia
        printf("%d -> ", rank);
        printf("0:");
        for (int j = 0; j < nr_proc1 - 1; j++) {
            printf("%d,", topologie[0][j]);
        }
        printf("%d 1:", topologie[0][nr_proc1 - 1]);
        for (int j = 0; j < nr_proc2 - 1; j++) {
            printf("%d,", topologie[1][j]);
        }
        printf("%d 2:", topologie[1][nr_proc2 - 1]);
        for (int j = 0; j < nr_proc3 - 1; j++) {
            printf("%d,", topologie[2][j]);
        }
        printf("%d", topologie[2][nr_proc3 - 1]);
        printf("\n");
    }   
    if (rank == 2) {
        FILE *clst2 = fopen("cluster2.txt", "r");
        fscanf(clst2, "%d", &nr_proc3);
        vec2 = malloc(nr_proc3 * sizeof(int));
        for (int i = 0; i < nr_proc3; i++) {
            fscanf(clst2, "%d", &vec2[i]);
            topologie[rank][i] = vec2[i];
        }
        // trimite nr de procese la coord 0
        MPI_Send(&nr_proc3, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        // trimite nr de procese la coord 1
        MPI_Send(&nr_proc3, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        // trimite workerii la coord 0
        MPI_Send(vec2, nr_proc3, MPI_INT, 0, 0, MPI_COMM_WORLD);
        // trimite workerii la coord 1
        MPI_Send(vec2, nr_proc3, MPI_INT, 1, 0, MPI_COMM_WORLD);

        // primeste nr de proc de la coord 0, aloca mem si primeste si workerii
        MPI_Recv(&nr_proc1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        vec0 = malloc(nr_proc1 * sizeof(int));
        MPI_Recv(vec0, nr_proc1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        for (int i = 0; i < nr_proc1; i++) {
            topologie[0][i] = vec0[i];
        }
        // primeste nr de proc de la coord 1, aloca mem si primeste si workerii
        MPI_Recv(&nr_proc2, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(1,%d)\n", rank);
        vec1 = malloc(nr_proc2 * sizeof(int));
        MPI_Recv(vec1, nr_proc2, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(1,%d)\n", rank);
        for (int i = 0; i < nr_proc2; i++) {
            topologie[1][i] = vec1[i];
        }
        // primesc de la 0 numarul si vectorul de calculat
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        calculateVector = malloc(N * sizeof(int));
        MPI_Recv(calculateVector, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);

        // primesc vectorul de start si de end 
        MPI_Recv(vecSCLS, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        MPI_Recv(vecECLS, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("M(0,%d)\n", rank);
        // prelucrez distanta pentru a calcula startul si endul fiecarui worker
        int aux = vecECLS[2] - vecSCLS[2];
        for (int i = 0; i < nr_proc3; i++) {
        vecStart[i] = (i *  aux / nr_proc3) + vecSCLS[2];
        vecEnd[i] = min(((i + 1) * aux) / nr_proc3, aux)  + vecSCLS[2];
        }

        // pentru toate procesele din coord 0 trimit toti workerii si vec de calculat
        for (int i = 0; i < nr_proc3; i++) {
             MPI_Send(&nr_proc1, 1, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
             MPI_Send(&nr_proc2, 1, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
             MPI_Send(&nr_proc3, 1, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);   
            for (int j = 0; j < 3; j++) {
                MPI_Send(topologie[j], procs - 3, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
            }
            MPI_Send(&vecStart[i], 1, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
            MPI_Send(&vecEnd[i], 1, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
            MPI_Send(calculateVector, N, MPI_INT, vec2[i], 0, MPI_COMM_WORLD);
        }
        // primesc vectorul schimbat si il modific in cel original pentru a-l trimite la coord 0
        int auxV[N];
        for (int i = 0; i < nr_proc3; i++) {
            MPI_Recv(auxV, N, MPI_INT, vec2[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("M(%d,%d)\n", vec2[i], rank);
            for (int j = vecStart[i]; j < vecEnd[i]; j++) {
                calculateVector[j] = auxV[j];
            }
        }
        // printf("2 -> Rezultat :");
        // for (int i = 0; i < N; i++) {
        //     printf("%d ", calculateVector[i]);
        // }
        // printf("\n");
        MPI_Send(calculateVector, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
        //printez topologia
        printf("%d -> ", rank);
        printf("0:");
        for (int j = 0; j < nr_proc1 - 1; j++) {
            printf("%d,", topologie[0][j]);
        }
        printf("%d 1:", topologie[0][nr_proc1 - 1]);
        for (int j = 0; j < nr_proc2 - 1; j++) {
            printf("%d,", topologie[1][j]);
        }
        printf("%d 2:", topologie[1][nr_proc2 - 1]);
        for (int j = 0; j < nr_proc3 - 1; j++) {
            printf("%d,", topologie[2][j]);
        }
        printf("%d", topologie[2][nr_proc3 - 1]);
        printf("\n");
        
    }
    int leader;
    if (rank != 0 && rank != 1 && rank != 2) {
        // primesc topologia
        MPI_Status status;
        MPI_Recv(&nr_proc1, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        leader = status.MPI_SOURCE;

        printf("M(%d,%d)\n", leader, rank);
        MPI_Recv(&nr_proc2, 1, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
        printf("M(%d,%d)\n", leader, rank);

        MPI_Recv(&nr_proc3, 1, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
        printf("M(%d,%d)\n", leader, rank);

        for (int j = 0; j < 3; j++) {
            MPI_Recv(topologie[j], procs - 3, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
            printf("M(%d,%d)\n", leader, rank);

        }
        // primesc vectorul de valori, il prelucrez si il trimit;
        MPI_Recv(&start, 1, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
        printf("M(%d,%d)\n", leader, rank);
        MPI_Recv(&end, 1, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
        printf("M(%d,%d)\n", leader, rank);
        MPI_Recv(&N, 1, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
        printf("M(%d,%d)\n", leader, rank);
        calculateVector = malloc(N * sizeof(int));
        MPI_Recv(calculateVector, N, MPI_INT, leader, 0, MPI_COMM_WORLD, &status);
        printf("M(%d,%d)\n", leader, rank);

        for (int i = start; i < end; i++) {
            calculateVector[i] *= 2;
        }
        MPI_Send(calculateVector, N, MPI_INT, leader, 0 , MPI_COMM_WORLD);

        // printez topologia pentru fiecare worker
        printf("%d -> ", rank);
        printf("0:");
        for (int j = 0; j < nr_proc1 - 1; j++) {
            printf("%d,", topologie[0][j]);
        }
        printf("%d 1:", topologie[0][nr_proc1 - 1]);
        for (int j = 0; j < nr_proc2 - 1; j++) {
            printf("%d,", topologie[1][j]);
        }
        printf("%d 2:", topologie[1][nr_proc2 - 1]);
        for (int j = 0; j < nr_proc3 - 1; j++) {
            printf("%d,", topologie[2][j]);
        }
        printf("%d", topologie[2][nr_proc3 - 1]);
        printf("\n");
    }
    MPI_Finalize();
}