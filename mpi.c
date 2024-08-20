#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define a structure to hold the coordinates of a point
struct Point {
    double x;
    double y;
};

// Function to read input file and populate points array
void read_input_file(const char *filename, double *D, int *MaxIterations, struct Point **points, int *num_points) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Cannot open file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read the first line to get D and MaxIterations values
    fscanf(file, "%lf %d", D, MaxIterations);

    int capacity = 9;
    *points = malloc(capacity * sizeof(struct Point));
    *num_points = 0;

    struct Point point;
    // Read the coordinates from the file and store them in the points array
    while (fscanf(file, "%lf %lf", &point.x, &point.y) == 2) {
        if (*num_points >= capacity) {
            capacity *= 2;
            *points = realloc(*points, capacity * sizeof(struct Point));
        }
        (*points)[(*num_points)++] = point;
    }

    fclose(file);
}

// Function to distribute points to all processes
void distribute_points(struct Point *points, struct Point *my_point, int rank) {
    // Scatter points to all processes
    MPI_Scatter(points, 2, MPI_DOUBLE, my_point, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("Rank %d received point: %lf, %lf\n", rank, my_point->x, my_point->y);
}

// Function to create a Cartesian topology for processes
MPI_Comm create_cartesian_topology(int K) {
    int dims[2] = {K, K};
    int periods[2] = {0, 0};  // No periodic neighbors
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    return cart_comm;
}

// Function to average the coordinates of neighbor points
struct Point average_points(struct Point *neighbors, int count) {
    double sum_x = 0, sum_y = 0;
    for (int i = 0; i < count; ++i) {
        sum_x += neighbors[i].x;
        sum_y += neighbors[i].y;
    }
    struct Point result = { sum_x / count, sum_y / count };
    return result;
}

// Function to calculate the distance of a point from the origin
double distance_to_origin(struct Point p) {
    return sqrt(p.x * p.x + p.y * p.y);
}

// Function to iterate over points and update their coordinates
void iterate_points(MPI_Comm cart_comm, struct Point *my_point, double D, int MaxIterations, int rank) {
    int up, down, left, right;
    // Determine the ranks of the neighboring processes
    MPI_Cart_shift(cart_comm, 0, 1, &left, &right);
    MPI_Cart_shift(cart_comm, 1, 1, &up, &down);

    for (int iter = 0; iter < MaxIterations; ++iter) {
        struct Point neighbors[4];
        int count = 0;

        // Communicate with left neighbor
        if (left != MPI_PROC_NULL) {
            MPI_Sendrecv(my_point, 2, MPI_DOUBLE, left, 0, &neighbors[count++], 2, MPI_DOUBLE, left, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        // Communicate with right neighbor
        if (right != MPI_PROC_NULL) {
            MPI_Sendrecv(my_point, 2, MPI_DOUBLE, right, 0, &neighbors[count++], 2, MPI_DOUBLE, right, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        // Communicate with up neighbor
        if (up != MPI_PROC_NULL) {
            MPI_Sendrecv(my_point, 2, MPI_DOUBLE, up, 0, &neighbors[count++], 2, MPI_DOUBLE, up, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        // Communicate with down neighbor
        if (down != MPI_PROC_NULL) {
            MPI_Sendrecv(my_point, 2, MPI_DOUBLE, down, 0, &neighbors[count++], 2, MPI_DOUBLE, down, 0, cart_comm, MPI_STATUS_IGNORE);
        }

        // Average the coordinates of the neighbors
        if (count > 0) {
            *my_point = average_points(neighbors, count);
        }

        double dist = distance_to_origin(*my_point);
        printf("Rank %d new point: %lf, %lf (distance: %lf)\n", rank, my_point->x, my_point->y, dist);

        // Check if all points are within the distance D
        int all_within_distance;
        int within_distance = (dist < D) ? 1 : 0;
        MPI_Allreduce(&within_distance, &all_within_distance, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

        // Terminate if all points are within the distance D
        if (all_within_distance) {
            if (rank == 0) {
                printf("Terminated after %d iterations\n", iter);
            }
            return;
        }
    }

    if (rank == 0) {
        printf("Stopped after maximum iterations\n");
    }
}

// Function to gather and print results
void output_results(struct Point *my_point, int rank, int size) {
    struct Point *all_points = NULL;
    if (rank == 0) {
        all_points = malloc(size * sizeof(struct Point));
    }

    MPI_Gather(my_point, 2, MPI_DOUBLE, all_points, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            printf("%lf %lf\n", all_points[i].x, all_points[i].y);
        }
        free(all_points);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double D;
    int MaxIterations;
    struct Point *points = NULL;
    int num_points;

    if (rank == 0) {
        // Read input file to get D, MaxIterations, and points array
        read_input_file("data.txt", &D, &MaxIterations, &points, &num_points);
        printf("D: %lf, MaxIterations: %d, NumPoints: %d\n", D, MaxIterations, num_points);
    }

    // Broadcast D and MaxIterations to all processes
    MPI_Bcast(&D, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&MaxIterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

    struct Point my_point;
    // Distribute points to all processes
    distribute_points(points, &my_point, rank);

    int K = (int)sqrt(size);
    // Create Cartesian topology
    MPI_Comm cart_comm = create_cartesian_topology(K);

    // Iterate over points to update their coordinates
    iterate_points(cart_comm, &my_point, D, MaxIterations, rank);

    // Gather and output the results
    output_results(&my_point, rank, size);

    if (rank == 0) {
        free(points);
    }

    MPI_Finalize();
    return 0;
}

