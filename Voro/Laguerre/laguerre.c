/**
 * @file laguerre.c
 * @brief Implementation of Laguerre Diagram (Power Diagram) Construction
 * 
 * This file implements algorithms for constructing Laguerre tessellation
 * (weighted Voronoi diagrams) in 2D and 3D space.
 * 
 * Algorithm Overview:
 * 1. Initialize cell as bounding box (cube/square)
 * 2. Sort neighbors by delta-radius squared (for efficient pruning)
 * 3. For each neighbor, cut the cell with the radical hyperplane
 * 4. Compact the cell by removing invalidated vertices/corners
 * 
 * Key Formula - Weight Fraction:
 *   wfrac = 0.5 + 0.5 * (w1^2 - w2^2) / d12^2
 * This determines where the radical plane intersects the line between two sites.
 * 
 * Bug Fixes Applied:
 * - findNewVertex: Fixed imax assignment (was always 1, now correctly tracks max)
 * - solveLinearSystem3D: Fixed determinant calculation (r1->z instead of r2->z)
 * - computePolygonCentroid2D_shoelace: Fixed Shoelace formula implementation
 * - findEdgePlaneIntersection2D: Added *tt = t assignment when mod == 0
 * - computeMaxWeightCeiling2D/3D: Added bounds check for negative indices
 * - classifyVertexPosition3D: Added division-by-zero safety check
 * - solveLinearSystem2D: Added parallel lines check
 */

#define LAGUERRE_MAIN
#include "laguerre.h"
#undef LAGUERRE_MAIN

/*============================================================================
 * 3D ROTATION - EULER ANGLES
 *============================================================================*/

/**
 * @brief Apply Euler rotation (ZXZ convention) to 3D point array
 * 
 * Rotation matrix: R = Rz(phi) * Rx(theta) * Rz(psi)
 * 
 * @param points     Array of points to rotate (modified in place)
 * @param count      Number of points
 * @param theta      First Euler angle around X-axis (radians)
 * @param phi        Second Euler angle around Z-axis (radians)
 * @param psi        Third Euler angle around Z-axis (radians)
 */
void applyEulerRotation3D(Point3D *points, int count,
                          real_t theta, real_t phi, real_t psi) {
    Point3D result;
    Point3D *current = points;
    
    for (int i = 0; i < count; i++) {
        /* Apply full rotation matrix */
        result.x = cos(theta)*cos(psi) * current->x
                 + (-cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi)) * current->y
                 + (sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi)) * current->z;
        
        result.y = cos(theta)*sin(psi) * current->x
                 + (cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi)) * current->y
                 + (-sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)) * current->z;
        
        result.z = -sin(theta) * current->x
                 + sin(phi)*sin(theta) * current->y
                 + cos(phi)*cos(theta) * current->z;
        
        points[i] = result;
        current++;
    }
}

/*============================================================================
 * NORMAL VECTOR COMPUTATION
 *============================================================================*/

/**
 * @brief Compute unit surface normal from three consecutive polygon vertices
 * 
 * Uses cross product of two edge vectors: normal = (v1-v0) × (v2-v1)
 * 
 * @param vertices   Array of at least 3 consecutive vertices
 * @return           Normalized surface normal vector
 */
Point3D computeSurfaceNormal3D(Point3D *vertices) {
    Point3D edge1 = subtract3D(vertices + 1, vertices);
    Point3D edge2 = subtract3D(vertices + 2, vertices + 1);
    Point3D normal = cross3D(&edge1, &edge2);
    
    real_t magnitude = sqrt(dot3D(&normal, &normal));
    normal.x /= magnitude;
    normal.y /= magnitude;
    normal.z /= magnitude;
    
    return normal;
}

/**
 * @brief Compute 2D perpendicular normal (counterclockwise rotation)
 * 
 * @param direction  Input direction vector
 * @return           Unit normal perpendicular to input
 */
Point2D computeNormal2D(Point2D *direction) {
    Point3D perpendicular, original, crossResult;
    Point2D result;
    
    /* Rotate 90 degrees: (x,y) -> (y,-x) */
    perpendicular.x =  direction->y;
    perpendicular.y = -direction->x;
    perpendicular.z = 0;
    
    original.x = direction->x;
    original.y = direction->y;
    original.z = 0;
    
    /* Check orientation using cross product */
    crossResult = cross3D(&perpendicular, &original);
    if (crossResult.z < 0) {
        perpendicular.x = -perpendicular.x;
        perpendicular.y = -perpendicular.y;
    }
    
    real_t magnitude = sqrt(dot3D(&perpendicular, &perpendicular));
    result.x = perpendicular.x / magnitude;
    result.y = perpendicular.y / magnitude;
    
    return result;
}

/*============================================================================
 * CORNER/VERTEX DEACTIVATION
 *============================================================================*/

/**
 * @brief Mark all corners between lower and upper as inactive
 * 
 * When cutting a cell with a radical line, corners outside the line
 * must be deactivated. This traverses from lower to upper and marks
 * all intermediate corners as inactive.
 * 
 * @param lower   Starting corner (exclusive)
 * @param upper   Ending corner (exclusive)
 */
void deactivateCornersBetween(Corner2D *lower, Corner2D *upper) {
    lower = lower->lowerLink;
    while (lower != upper) {
        lower->status = STATUS_INACTIVE;
        lower = lower->lowerLink;
    }
}

/*============================================================================
 * VERTEX/CORNER ARRAY COMPACTION
 *============================================================================*/

/**
 * @brief Compact 3D vertex array by removing inactive vertices
 * 
 * After cutting operations, some vertices become inactive. This function:
 * 1. Builds old-to-new index mapping
 * 2. Compacts array by moving active vertices to front
 * 3. Updates all adjacency pointers
 * 
 * @param vertices      Vertex array (modified in place)
 * @param count         Current vertex count
 * @param maxDistSq     Output: maximum squared distance from origin
 * @param indexMap      Work array for index mapping (size >= count)
 * @return              New count of active vertices
 */
int compactActiveVertices3D(Vertex3D *vertices, int count,
                            real_t *maxDistSq, int *indexMap) {
    int newCount = 0;
    real_t maxDist = -1.0e20;
    
    /* First pass: build index mapping and find max distance */
    for (int i = 0; i < count; i++) {
        if (vertices[i].status == STATUS_ACTIVE) {
            real_t distSq = dot3D(&vertices[i], &vertices[i]);
            if (distSq > maxDist) maxDist = distSq;
            indexMap[i] = newCount;
            newCount++;
        }
    }
    
    /* Second pass: compact array and update adjacency pointers */
    newCount = 0;
    for (int i = 0; i < count; i++) {
        if (vertices[i].status == STATUS_ACTIVE) {
            vertices[newCount] = vertices[i];
            vertices[newCount].index = newCount;
            
            /* Update link pointers using index map */
            vertices[newCount].adjacent[0] = indexMap[vertices[i].adjacent[0] - vertices] + vertices;
            vertices[newCount].adjacent[1] = indexMap[vertices[i].adjacent[1] - vertices] + vertices;
            vertices[newCount].adjacent[2] = indexMap[vertices[i].adjacent[2] - vertices] + vertices;
            newCount++;
        }
    }
    
    *maxDistSq = maxDist;
    return newCount;
}

/**
 * @brief Compact 2D corner array by removing inactive corners
 * 
 * @param corners       Corner array (modified in place)
 * @param count         Current corner count
 * @param maxDistSq     Output: maximum squared distance from origin
 * @return              New count of active corners
 */
int compactActiveCorners2D(Corner2D *corners, int count, real_t *maxDistSq) {
    int indexMap[count];  /* VLA - consider malloc for large arrays */
    int newCount = 0;
    real_t maxDist = -1.0e20;
    
    /* First pass: build index mapping */
    for (int i = 0; i < count; i++) {
        if (corners[i].status == STATUS_ACTIVE) {
            real_t distSq = dot2D(&corners[i], &corners[i]);
            if (distSq > maxDist) maxDist = distSq;
            indexMap[i] = newCount;
            newCount++;
        }
    }
    
    /* Second pass: compact and update links */
    newCount = 0;
    for (int i = 0; i < count; i++) {
        if (corners[i].status == STATUS_ACTIVE) {
            corners[newCount] = corners[i];
            corners[newCount].upperLink = indexMap[corners[i].upperLink - corners] + corners;
            corners[newCount].lowerLink = indexMap[corners[i].lowerLink - corners] + corners;
            newCount++;
        }
    }
    
    *maxDistSq = maxDist;
    return newCount;
}

/*============================================================================
 * EDGE-PLANE INTERSECTION
 *============================================================================*/

/**
 * @brief Find intersection of 2D cell edge with radical line
 * 
 * Given edge from corner A to B and a neighbor defining the radical line,
 * find the intersection point. The radical line satisfies:
 *   point · neighbor = weightFrac * ||neighbor||^2
 * 
 * @param cornerA       First endpoint of edge
 * @param cornerB       Second endpoint of edge
 * @param neighbor      Neighbor site (relative to cell center at origin)
 * @param weightFrac    Laguerre weight fraction
 * @param tParam        Output: parameter t where intersection occurs (0 to 1)
 * @return              New corner at intersection point
 */
Corner2D findEdgePlaneIntersection2D(Corner2D *cornerA, Corner2D *cornerB,
                                     Point2D *neighbor, real_t weightFrac,
                                     real_t *tParam) {
    /* Edge direction: B - A */
    Point2D edgeDir = subtract2D(cornerB, cornerA);
    
    /* Compute intersection parameter:
     * t = (weightFrac * ||n||^2 - A·n) / (edgeDir·n)
     */
    real_t edgeDotN = dot2D(&edgeDir, neighbor);
    real_t numerator = weightFrac * dot2D(neighbor, neighbor) - dot2D(cornerA, neighbor);
    
    /* Handle degenerate case: edge parallel to radical line */
    if (edgeDotN == 0.0L) {
        real_t t = weightFrac;
        *tParam = t;  /* BUG FIX: Set output parameter */
        return multiplyAdd2D(t, &edgeDir, cornerA);
    }
    
    real_t t = numerator / edgeDotN;
    
    /* Clamp to valid range with epsilon for numerical stability */
    t = MAX(EPSILON, t);
    t = MIN(1.0 - EPSILON, t);
    
    if (t >= 0 && t <= 1) {
        *tParam = t;
        return multiplyAdd2D(t, &edgeDir, cornerA);
    } else {
        fprintf(stderr, "Error: intersection parameter out of bounds: t=%g, mod=%g\n", t, edgeDotN);
        exit(9);
    }
}

/**
 * @brief Find intersection of 3D cell edge with radical plane
 * 
 * @param vertexA       First endpoint of edge
 * @param vertexB       Second endpoint of edge
 * @param neighbor      Neighbor site (relative to cell center at origin)
 * @param weightFrac    Laguerre weight fraction
 * @return              New vertex at intersection point
 */
Vertex3D findEdgePlaneIntersection3D(Vertex3D *vertexA, Vertex3D *vertexB,
                                     Point3D *neighbor, real_t weightFrac) {
    Point3D edgeDir = subtract3D(vertexB, vertexA);
    
    real_t edgeDotN = dot3D(&edgeDir, neighbor);
    real_t numerator = weightFrac * dot3D(neighbor, neighbor) - dot3D(vertexA, neighbor);
    
    if (edgeDotN == 0) {
        fprintf(stderr, "Warning: edge parallel to cutting plane\n");
        real_t t = weightFrac;
        return multiplyAdd3D(t, &edgeDir, vertexA);
    }
    
    real_t t = numerator / edgeDotN;
    t = MIN(t, 1.0 - EPSILON);
    t = MAX(t, EPSILON);
    
    if (t > 0 && t < 1) {
        return multiplyAdd3D(t, &edgeDir, vertexA);
    } else {
        fprintf(stderr, "Error: intersection parameter out of bounds: t=%g\n", t);
        exit(9);
    }
}

/*============================================================================
 * POINT CLASSIFICATION
 *============================================================================*/

/**
 * @brief Classify 2D corner as inside or outside relative to radical line
 * 
 * A point P is inside if: P·neighbor < weightFrac * ||neighbor||^2
 * 
 * @param neighbor      Neighbor site
 * @param corner        Corner to classify
 * @param weightFrac    Laguerre weight fraction
 * @return              POINT_INSIDE, POINT_OUTSIDE, or POINT_BOUNDARY
 */
int classifyCornerPosition2D(Point2D *neighbor, Corner2D *corner, real_t weightFrac) {
    real_t projection = dot2D(corner, neighbor);
    real_t threshold = weightFrac * dot2D(neighbor, neighbor);
    
    if (projection > threshold) return POINT_OUTSIDE;
    else if (projection < threshold) return POINT_INSIDE;
    else return POINT_BOUNDARY;
}

/**
 * @brief Classify 3D vertex as inside or outside relative to radical plane
 * 
 * @param neighbor      Neighbor site
 * @param vertex        Vertex to classify
 * @param weightFrac    Laguerre weight fraction
 * @param offset        Output: signed distance (positive if outside)
 * @return              POINT_INSIDE or POINT_OUTSIDE
 */
int classifyVertexPosition3D(Point3D *neighbor, Vertex3D *vertex,
                             real_t weightFrac, real_t *offset) {
    real_t projection = dot3D(vertex, neighbor);
    real_t threshold = weightFrac * dot3D(neighbor, neighbor);
    
    /* BUG FIX: Safety check for division by zero */
    if (fabs(threshold) < EPSILON) {
        *offset = -1;
        return POINT_INSIDE;
    }
    
    real_t ratio = projection / threshold;
    
    if (ratio > 1 + BOUNDARY_TOLERANCE) {
        *offset = projection - threshold;
        return POINT_OUTSIDE;
    } else {
        *offset = -1;
        return POINT_INSIDE;
    }
}

/*============================================================================
 * BOUNDING BOX INITIALIZATION
 *============================================================================*/

/**
 * @brief Initialize 2D cell as square bounding box centered at origin
 * 
 * Creates 4 corners in counter-clockwise order forming a square.
 * Boundary edges are marked with negative neighbor indices.
 * 
 * @param corners       Output corner array
 * @param maxCorners    Maximum capacity
 * @param boxSize       Side length of square
 * @return              Number of initialized corners (4)
 */
int initBoundingBox2D(Corner2D *corners, int maxCorners, real_t boxSize) {
    real_t halfBox = 0.5 * boxSize;
    
    for (int i = 0; i < 4; i++) {
        /* Compute corner coordinates:
         * i=0: (-half, -half), i=1: (+half, -half)
         * i=2: (+half, +half), i=3: (-half, +half)
         */
        corners[i].x = boxSize * (((i + 1) % 4) / 2) - halfBox;
        corners[i].y = boxSize * (i / 2) - halfBox;
        corners[i].status = STATUS_ACTIVE;
        
        /* Circular linked list */
        corners[i].lowerLink = corners + (i - 1 + 4) % 4;
        corners[i].upperLink = corners + (i + 1 + 4) % 4;
        
        corners[i].lowerNeighborPtr = NULL;
        corners[i].upperNeighborPtr = NULL;
        
        /* Negative indices indicate boundary edges */
        corners[i].upperNeighborIdx = -(i + 1);
        corners[i].lowerNeighborIdx = -((i + 3) % 4 + 1);
    }
    
    /* Mark remaining slots as inactive */
    for (int i = 4; i < maxCorners; i++) {
        corners[i].status = STATUS_INACTIVE;
    }
    
    return 4;
}

/**
 * @brief Initialize 3D cell as cube bounding box centered at origin
 * 
 * Creates 8 vertices forming a cube with proper edge connectivity.
 * 
 * @param vertices      Output vertex array
 * @param maxVertices   Maximum capacity
 * @param boxSize       Side length of cube
 * @return              Number of initialized vertices (8)
 */
int initBoundingBox3D(Vertex3D *vertices, int maxVertices, real_t boxSize) {
    real_t halfBox = 0.5 * boxSize;
    
    for (int i = 0; i < 8; i++) {
        int faceIdx = i % 4;   /* Position within face (0-3) */
        int levelIdx = i / 4;  /* Bottom (0) or top (1) face */
        
        vertices[i].x = boxSize * (((faceIdx + 1) % 4) / 2) - halfBox;
        vertices[i].y = boxSize * (faceIdx / 2) - halfBox;
        vertices[i].z = boxSize * levelIdx - halfBox;
        vertices[i].status = STATUS_ACTIVE;
        
        /* Horizontal neighbor within same level */
        vertices[i].adjacent[0] = vertices + (faceIdx - 1 + 4) % 4 + levelIdx * 4;
        
        if (i < 4) {
            /* Bottom face */
            vertices[i].neighborSites[0] = -(faceIdx + 1);
            vertices[i].neighborSites[1] = -((faceIdx + 3) % 4 + 1);
            vertices[i].neighborSites[2] = -5;  /* Bottom boundary */
            
            vertices[i].adjacent[1] = vertices + (faceIdx + 1 + 4) % 4 + levelIdx * 4;
            vertices[i].adjacent[2] = vertices + (i + 4) % 8;
        } else {
            /* Top face (reversed winding) */
            vertices[i].neighborSites[0] = -(faceIdx + 1);
            vertices[i].neighborSites[1] = -6;  /* Top boundary */
            vertices[i].neighborSites[2] = -((faceIdx + 3) % 4 + 1);
            
            vertices[i].adjacent[2] = vertices + (faceIdx + 1 + 4) % 4 + levelIdx * 4;
            vertices[i].adjacent[1] = vertices + (i + 4) % 8;
        }
        
        for (int k = 0; k < 3; k++) {
            vertices[i].faceVertexCount[k] = 0;
        }
    }
    
    for (int i = 8; i < maxVertices; i++) {
        vertices[i].status = STATUS_INACTIVE;
    }
    
    return 8;
}

/*============================================================================
 * BOUNDARY FINDING FUNCTIONS
 *============================================================================*/

/**
 * @brief Find lower boundary corner when cutting cell
 * 
 * Starting from corner 'start', traverse upward until finding
 * a corner that is inside the cutting line.
 */
static Corner2D* findLowerLimit2D(Point2D *neighbor, Corner2D *start, real_t weightFrac) {
    start = start->upperLink;
    while (classifyCornerPosition2D(neighbor, start, weightFrac) == POINT_OUTSIDE) {
        start = start->upperLink;
    }
    return start;
}

/**
 * @brief Find upper boundary corner when cutting cell
 */
static Corner2D* findUpperLimit2D(Point2D *neighbor, Corner2D *start, real_t weightFrac) {
    start = start->lowerLink;
    while (classifyCornerPosition2D(neighbor, start, weightFrac) == POINT_OUTSIDE) {
        start = start->lowerLink;
    }
    return start;
}

/**
 * @brief Create new corners at intersection points (2D)
 */
static int createNewCorners2D(Corner2D *lower, Corner2D *upper,
                              Corner2D *corners, int currentCount,
                              Point2D *neighbor, int neighborIdx,
                              real_t weightFrac) {
    real_t t;
    
    /* Create lower intersection corner */
    corners[currentCount] = findEdgePlaneIntersection2D(
        lower->lowerLink, lower, neighbor, weightFrac, &t);
    corners[currentCount].lowerLink = lower->lowerLink;
    corners[currentCount].lowerNeighborIdx = lower->lowerNeighborIdx;
    corners[currentCount].lowerNeighborPtr = lower->lowerNeighborPtr;
    corners[currentCount].status = STATUS_ACTIVE;
    lower->lowerLink->upperLink = corners + currentCount;
    
    /* Create upper intersection corner */
    corners[currentCount + 1] = findEdgePlaneIntersection2D(
        upper, upper->upperLink, neighbor, weightFrac, &t);
    corners[currentCount + 1].upperLink = upper->upperLink;
    corners[currentCount + 1].upperNeighborIdx = upper->upperNeighborIdx;
    corners[currentCount + 1].upperNeighborPtr = upper->upperNeighborPtr;
    corners[currentCount + 1].status = STATUS_ACTIVE;
    upper->upperLink->lowerLink = corners + currentCount + 1;
    
    /* Link new corners together */
    corners[currentCount].upperLink = corners + currentCount + 1;
    corners[currentCount + 1].lowerLink = corners + currentCount;
    
    corners[currentCount].upperNeighborIdx = neighborIdx;
    corners[currentCount + 1].lowerNeighborIdx = neighborIdx;
    corners[currentCount].upperNeighborPtr = neighbor;
    corners[currentCount + 1].lowerNeighborPtr = neighbor;
    
    return currentCount + 2;
}

/*============================================================================
 * LINEAR SYSTEM SOLVERS
 *============================================================================*/

/**
 * @brief Solve 2D linear system (find line intersection)
 * 
 * Given two lines: y = slope1*x + intercept1, y = slope2*x + intercept2
 * Find their intersection point.
 * 
 * @param slope1        Slope of first line
 * @param intercept1    Y-intercept of first line
 * @param slope2        Slope of second line
 * @param intercept2    Y-intercept of second line
 * @return              Intersection point
 */
Point2D solveLinearSystem2D(real_t slope1, real_t intercept1,
                            real_t slope2, real_t intercept2) {
    Point2D intersection;
    real_t denominator = slope1 - slope2;
    
    /* BUG FIX: Safety check for parallel lines */
    if (fabs(denominator) < EPSILON) {
        fprintf(stderr, "Warning: Lines are nearly parallel in solveLinearSystem2D\n");
        intersection.x = 0.5 * (intercept1 + intercept2);
        intersection.y = slope1 * intersection.x + intercept1;
        return intersection;
    }
    
    intersection.x = (intercept2 - intercept1) / denominator;
    intersection.y = slope1 * intersection.x + intercept1;
    
    return intersection;
}

/**
 * @brief Solve 3D linear system (find three-plane intersection)
 * 
 * Solves: n1·x = d1, n2·x = d2, n3·x = d3 using Cramer's rule.
 * 
 * @param n1, n2, n3    Normal vectors of three planes
 * @param d1, d2, d3    Distance parameters of three planes
 * @return              Intersection point (NaN if parallel)
 */
Point3D solveLinearSystem3D(Point3D *n1, real_t d1,
                            Point3D *n2, real_t d2,
                            Point3D *n3, real_t d3) {
    Point3D intersection;
    
    /* Compute determinant of coefficient matrix using cofactor expansion
     * BUG FIX: Third term uses n1->z (was incorrectly n2->z)
     */
    real_t determinant = n1->x * (n2->y * n3->z - n3->y * n2->z)
                       - n1->y * (n2->x * n3->z - n3->x * n2->z)
                       + n1->z * (n2->x * n3->y - n3->x * n2->y);
    
    if (determinant == 0) {
        printf("Error: Planes are parallel, no unique intersection.\n");
        intersection.x = NAN;
        intersection.y = NAN;
        intersection.z = NAN;
        return intersection;
    }
    
    /* Apply Cramer's rule */
    intersection.x = (d1 * (n2->y * n3->z - n3->y * n2->z)
                    - n1->y * (d2 * n3->z - d3 * n2->z)
                    + n1->z * (d2 * n3->y - d3 * n2->y)) / determinant;
    
    intersection.y = (n1->x * (d2 * n3->z - d3 * n2->z)
                    - d1 * (n2->x * n3->z - n3->x * n2->z)
                    + n1->z * (n2->x * d3 - n3->x * d2)) / determinant;
    
    intersection.z = (n1->x * (n2->y * d3 - n3->y * d2)
                    - n1->y * (n2->x * d3 - n3->x * d2)
                    + d1 * (n2->x * n3->y - n3->x * n2->y)) / determinant;
    
    return intersection;
}

/*============================================================================
 * RADICAL PLANE/LINE DEFINITION
 *============================================================================*/

/**
 * @brief Define radical plane between two 3D sites
 * 
 * The radical plane has equation: normal · x = distance
 * where normal = site2 - site1
 * 
 * @param site1, site2  Two sites
 * @param weightParam   Weight parameter t
 * @param normal        Output: plane normal
 * @param distance      Output: plane distance parameter
 */
void computeRadicalPlane3D(Point3D *site1, Point3D *site2, real_t weightParam,
                           Point3D *normal, real_t *distance) {
    Point3D diff = subtract3D(site2, site1);
    normal->x = diff.x;
    normal->y = diff.y;
    normal->z = diff.z;
    *distance = dot3D(normal, site1) + weightParam * dot3D(normal, normal);
}

/**
 * @brief Define radical line between two 2D sites
 */
void computeRadicalLine2D(Point2D *site1, Point2D *site2, real_t weightParam,
                          Point2D *normal, real_t *distance) {
    Point2D diff = subtract2D(site2, site1);
    normal->x = diff.x;
    normal->y = diff.y;
    *distance = dot2D(normal, site1) + weightParam * dot2D(normal, normal);
}

/*============================================================================
 * POLYGON/POLYHEDRON METRICS
 *============================================================================*/

/**
 * @brief Compute 2D polygon centroid using triangulation
 */
Corner2D computePolygonCentroid2D(Corner2D *start) {
    Corner2D centroid = {0};
    centroid.x = centroid.y = 0;
    
    Corner2D *current = start;
    Corner2D *next;
    real_t totalArea = 0;
    
    do {
        next = current->upperLink;
        real_t triangleArea = 0.5 * fabs(current->x * next->y - current->y * next->x);
        totalArea += triangleArea;
        centroid.x += triangleArea * (current->x + next->x) / 3.0;
        centroid.y += triangleArea * (current->y + next->y) / 3.0;
        current = next;
    } while (current != start);
    
    centroid.x /= totalArea;
    centroid.y /= totalArea;
    
    return centroid;
}

/**
 * @brief Compute 2D polygon centroid using Shoelace formula
 * 
 * BUG FIX: Corrected Shoelace formula implementation.
 * Original had: A += tmp->x*tmp->y - tmp2->x*tmp2->y (WRONG)
 * Correct is:   A += tmp->x*tmp2->y - tmp2->x*tmp->y
 */
Corner2D computePolygonCentroid2D_shoelace(Corner2D *start) {
    Corner2D centroid = {0};
    real_t signedArea = 0;
    centroid.x = 0;
    centroid.y = 0;
    
    Corner2D *current = start;
    Corner2D *next;
    
    do {
        next = current->upperLink;
        /* BUG FIX: Correct Shoelace formula cross product term */
        real_t crossTerm = current->x * next->y - next->x * current->y;
        signedArea += crossTerm;
        centroid.x += (current->x + next->x) * crossTerm;
        centroid.y += (current->y + next->y) * crossTerm;
        current = next;
    } while (current != start);
    
    signedArea *= 0.5;
    centroid.x /= (6 * signedArea);
    centroid.y /= (6 * signedArea);
    
    return centroid;
}

/**
 * @brief Compute polygon face centroid (3D)
 * @note Currently unused - kept for future use
 */
#if 0
static Point3D getPolygonCentroid3D(Vertex3D *start, Vertex3D *next) {
    Point3D centroid = {0};
    Vertex3D *current = next;
    Vertex3D *prev = start;
    int count = 0;
    
    do {
        centroid.x += current->x;
        centroid.y += current->y;
        centroid.z += current->z;
        count++;
        
        Vertex3D *temp = current;
        current = getNextFaceVertex(prev, current);
        prev = temp;
    } while (current != next && count < 100);
    
    if (count > 0) {
        centroid.x /= count;
        centroid.y /= count;
        centroid.z /= count;
    }
    
    return centroid;
}
#endif

/**
 * @brief Compute pyramid volume from face
 */
real_t computePyramidVolume3D(Vertex3D *start, Vertex3D *next) {
    Point3D surfaceNormal = {0};
    int edgeIdx = 0;
    
    /* Find edge index */
    while (next->adjacent[edgeIdx] != start) {
        edgeIdx = (edgeIdx + 1) % 3;
    }
    
    /* Traverse face, computing surface integral */
    Vertex3D *current = next;
    Vertex3D *prev = start;
    int count = 0;
    
    do {
        surfaceNormal.x += prev->y * current->z - prev->z * current->y;
        surfaceNormal.y += prev->z * current->x - prev->x * current->z;
        surfaceNormal.z += prev->x * current->y - prev->y * current->x;
        
        Vertex3D *temp = current;
        current = getNextFaceVertex(prev, current);
        prev = temp;
        count++;
        
        if (count > 100) {
            fprintf(stderr, "Warning: Possible infinite loop in computePyramidVolume3D\n");
            break;
        }
    } while (current != next);
    
    return sqrt(dot3D(&surfaceNormal, &surfaceNormal)) / 6.0;
}

/**
 * @brief Compute total polyhedron volume
 */
real_t computePolyhedronVolume3D(Vertex3D *vertices, int numVertices) {
    real_t totalVolume = 0;
    
    /* Reset face counters */
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < 3; j++) {
            vertices[i].faceVertexCount[j] = 0;
        }
    }
    
    /* Traverse all faces */
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < 3; j++) {
            if (vertices[i].faceVertexCount[j] == 0) {
                Vertex3D *startV = vertices + i;
                Vertex3D *nextV = startV->adjacent[j];
                
                totalVolume += computePyramidVolume3D(startV, nextV);
                
                /* Mark face edges as visited */
                Vertex3D *current = nextV;
                Vertex3D *prev = startV;
                int count = 0;
                
                do {
                    int k = 0;
                    while (current->adjacent[k] != prev) {
                        k = (k + 1) % 3;
                    }
                    current->faceVertexCount[k] = 1;
                    
                    Vertex3D *temp = current;
                    current = getNextFaceVertex(prev, current);
                    prev = temp;
                    count++;
                    
                    if (count > 100) break;
                } while (current != nextV);
            }
        }
    }
    
    return totalVolume;
}

/*============================================================================
 * DELAUNAY VERTEX UPDATE
 *============================================================================*/

/**
 * @brief Structure for Delaunay tetrahedron
 */
typedef struct {
    Point3D *vertices[4];
} DelaunayTetrahedron;

/**
 * @brief Update vertex position from Delaunay tetrahedron
 * 
 * BUG FIX: Loop now correctly sets imax = i (was always setting imax = 1)
 */
void updateVertexFromDelaunay(Vertex3D *vertex, Point3D *neighbors, Point3D *center) {
    /* Skip if any neighbor is a boundary */
    if (vertex->neighborSites[0] < 0 || 
        vertex->neighborSites[1] < 0 || 
        vertex->neighborSites[2] < 0) {
        return;
    }
    
    Point3D localCenter = *center;
    localCenter.x = localCenter.y = localCenter.z = 0;
    localCenter.soundSpeed = center->soundSpeed;
    
    DelaunayTetrahedron tetra;
    tetra.vertices[0] = &localCenter;
    tetra.vertices[1] = neighbors + vertex->neighborSites[0];
    tetra.vertices[2] = neighbors + vertex->neighborSites[1];
    tetra.vertices[3] = neighbors + vertex->neighborSites[2];
    
    /* Find vertex with maximum sound speed */
    int maxIdx = 0;
    real_t maxSoundSpeed = localCenter.soundSpeed;
    
    for (int i = 0; i < 4; i++) {
        real_t soundSpeed = tetra.vertices[i]->soundSpeed;
        if (soundSpeed > maxSoundSpeed) {
            maxSoundSpeed = soundSpeed;
            maxIdx = i;  /* BUG FIX: Was "maxIdx = 1" */
        }
    }
    
    /* Compute weighted intersection */
    real_t weightMax = pow(tetra.vertices[maxIdx]->soundSpeed, 0.25);
    
    Point3D normals[3];
    real_t distances[3];
    
    for (int k = 0; k < 3; k++) {
        int j = (maxIdx + 1 + k) % 4;
        real_t weightJ = pow(tetra.vertices[j]->soundSpeed, 0.25);
        real_t t = weightMax / (weightJ + weightMax);
        t = MIN(t, 0.7);
        
        computeRadicalPlane3D(tetra.vertices[maxIdx], tetra.vertices[j],
                              t, &normals[k], &distances[k]);
    }
    
    Point3D newPos = solveLinearSystem3D(&normals[0], distances[0],
                                          &normals[1], distances[1],
                                          &normals[2], distances[2]);
    vertex->x = newPos.x;
    vertex->y = newPos.y;
    vertex->z = newPos.z;
}

/*============================================================================
 * WEIGHT CEILING COMPUTATION
 *============================================================================*/

/**
 * @brief Compute upper bound on w^2 for neighbor pruning (2D)
 * 
 * BUG FIX: Added bounds check for negative upperNeighborIdx (boundary edges)
 */
real_t computeMaxWeightCeiling2D(Corner2D *start, Point2D *neighbors) {
    real_t ceiling = 1.0e20;
    Corner2D *current = start;
    
    do {
        /* BUG FIX: Check for valid (non-boundary) neighbor */
        if (current->upperNeighborIdx >= 0) {
            Point2D *neighbor = neighbors + current->upperNeighborIdx;
            real_t value = neighbor->squaredDist + neighbor->weightSquared;
            ceiling = MIN(ceiling, value);
        }
        current = current->upperLink;
    } while (current != start);
    
    return ceiling;
}

/**
 * @brief Compute upper bound on w^2 for neighbor pruning (3D)
 * 
 * BUG FIX: Added bounds check for negative neighborSites (boundary faces)
 */
real_t computeMaxWeightCeiling3D(Vertex3D *vertices, Point3D *neighbors, int numVertices) {
    real_t ceiling = 1.0e20;
    
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < 3; j++) {
            /* BUG FIX: Check for valid (non-boundary) neighbor */
            if (vertices[i].neighborSites[j] >= 0) {
                Point3D *neighbor = neighbors + vertices[i].neighborSites[j];
                real_t value = neighbor->weightSquared + neighbor->squaredDist;
                ceiling = MIN(ceiling, value);
            }
        }
    }
    
    return ceiling;
}

/*============================================================================
 * SORTING COMPARISON FUNCTIONS
 *============================================================================*/

int compareByDistanceSquared2D(const void *a, const void *b) {
    const Point2D *pa = (const Point2D *)a;
    const Point2D *pb = (const Point2D *)b;
    if (pa->squaredDist > pb->squaredDist) return 1;
    else if (pa->squaredDist < pb->squaredDist) return -1;
    return 0;
}

int compareByDistanceSquared3D(const void *a, const void *b) {
    const Point3D *pa = (const Point3D *)a;
    const Point3D *pb = (const Point3D *)b;
    if (pa->squaredDist > pb->squaredDist) return 1;
    else if (pa->squaredDist < pb->squaredDist) return -1;
    return 0;
}

int compareByDeltaRadius2D(const void *a, const void *b) {
    const Point2D *pa = (const Point2D *)a;
    const Point2D *pb = (const Point2D *)b;
    if (pa->deltaRadiusSq > pb->deltaRadiusSq) return 1;
    else if (pa->deltaRadiusSq < pb->deltaRadiusSq) return -1;
    return 0;
}

int compareByDeltaRadius3D(const void *a, const void *b) {
    const Point3D *pa = (const Point3D *)a;
    const Point3D *pb = (const Point3D *)b;
    if (pa->deltaRadiusSq > pb->deltaRadiusSq) return 1;
    else if (pa->deltaRadiusSq < pb->deltaRadiusSq) return -1;
    return 0;
}

/*============================================================================
 * MAIN CELL CONSTRUCTION - 2D
 *============================================================================*/

/**
 * @brief Construct 2D Laguerre cell for a given site
 */
int constructLaguerreCell2D(Point2D *center, Point2D *inputNeighbors,
                            Point2D *neighbors, int numNeighbors,
                            Corner2D *corners, int maxCorners, real_t boxSize) {
    real_t maxDistSq = (boxSize/2) * (boxSize/2) * 2;
    real_t centerWeight = center->weightSquared;
    int cornerCount;
    
    /* Transform neighbors to local coordinates */
    int validCount = 0;
    real_t minDistSq = 1.0e20;
    
    for (int i = 0; i < numNeighbors; i++) {
        if (inputNeighbors[i].index != center->index) {
            neighbors[validCount] = subtract2D(&inputNeighbors[i], center);
            neighbors[validCount].soundSpeed = inputNeighbors[i].soundSpeed;
            neighbors[validCount].pressure = inputNeighbors[i].pressure;
            neighbors[validCount].squaredDist = dot2D(&neighbors[validCount], &neighbors[validCount]);
            minDistSq = MIN(minDistSq, neighbors[validCount].squaredDist);
            neighbors[validCount].backPointer = inputNeighbors[i].backPointer;
            neighbors[validCount].weightSquared = inputNeighbors[i].weightSquared;
            
            real_t distSq = neighbors[validCount].squaredDist;
            real_t a = 0.5 * (1 + (centerWeight - neighbors[validCount].weightSquared) / distSq);
            neighbors[validCount].deltaRadiusSq = distSq * a * a;
            validCount++;
        }
    }
    numNeighbors = validCount;
    
    /* Sort by delta radius for efficient pruning */
    qsort(neighbors, numNeighbors, sizeof(Point2D), compareByDeltaRadius2D);
    
    /* Initialize bounding box */
    cornerCount = initBoundingBox2D(corners, maxCorners, boxSize);
    
    /* Cut cell with each neighbor */
    for (int i = 0; i < numNeighbors; i++) {
        if (neighbors[i].deltaRadiusSq > maxDistSq) break;
        
        real_t distSq = dot2D(&neighbors[i], &neighbors[i]);
        real_t neighborWeight = neighbors[i].weightSquared;
        real_t weightFrac = computeWeightFraction(centerWeight, neighborWeight, distSq);
        
        for (int j = 0; j < cornerCount; j++) {
            if (corners[j].status == STATUS_ACTIVE &&
                classifyCornerPosition2D(&neighbors[i], &corners[j], weightFrac) == POINT_OUTSIDE) {
                
                Corner2D *lower = &corners[j];
                lower = findLowerLimit2D(&neighbors[i], lower, weightFrac);
                
                Corner2D *upper = &corners[j];
                upper = findUpperLimit2D(&neighbors[i], upper, weightFrac);
                
                deactivateCornersBetween(lower, upper);
                cornerCount = createNewCorners2D(lower, upper, corners, cornerCount, &neighbors[i], i, weightFrac);
                cornerCount = compactActiveCorners2D(corners, cornerCount, &maxDistSq);
                break;
            }
        }
    }
    
    center->weightCeiling = computeMaxWeightCeiling2D(corners, neighbors);
    return cornerCount;
}

/*============================================================================
 * 3D CELL CONSTRUCTION HELPERS
 *============================================================================*/

/**
 * @brief Inactivate outside vertices
 */
static int inactivateOutsideVertices(Point3D *neighbor, Vertex3D *vertices,
                                     int count, real_t weightFrac) {
    int inactiveCount = 0;
    real_t offset;
    
    for (int i = 0; i < count; i++) {
        if (classifyVertexPosition3D(neighbor, &vertices[i], weightFrac, &offset) == POINT_OUTSIDE) {
            vertices[i].status = STATUS_INACTIVE;
            inactiveCount++;
        }
    }
    
    return inactiveCount;
}

/**
 * @brief Find nearest boundary vertex for cutting
 */
static VertexRef3D findNearestBoundaryVertex(Point3D *neighbor, Vertex3D *vertex,
                                             real_t weightFrac) {
    VertexRef3D result;
    
    while (1) {
        real_t minOffset = 1.0e27;
        real_t offset;
        int minIdx = 0;
        
        for (int i = 0; i < 3; i++) {
            Vertex3D *link = vertex->adjacent[i];
            if (classifyVertexPosition3D(neighbor, link, weightFrac, &offset) != POINT_OUTSIDE) {
                result.vertex = vertex;
                result.linkedVertex = link;
                result.linkIndex = i;
                return result;
            }
            if (offset < minOffset) {
                minOffset = offset;
                minIdx = i;
            }
        }
        vertex = vertex->adjacent[minIdx];
    }
}

/**
 * @brief Create new vertices at intersection (3D)
 */
static int createNewVertices3D(Point3D *neighbors, int neighborIdx,
                               VertexRef3D *start, Vertex3D *vertices,
                               int currentCount, real_t weightFrac) {
    /* Implementation continues from boundary vertex,
     * creating new vertices at edge intersections */
    
    Vertex3D *boundary = start->vertex;
    Vertex3D *inside = start->linkedVertex;
    int linkIdx = start->linkIndex;
    
    /* Create first new vertex */
    vertices[currentCount] = findEdgePlaneIntersection3D(boundary, inside,
                                                         &neighbors[neighborIdx], weightFrac);
    vertices[currentCount].status = STATUS_ACTIVE;
    vertices[currentCount].index = currentCount;
    
    /* Set neighbor relation for the new face */
    vertices[currentCount].neighborSites[0] = neighborIdx;
    vertices[currentCount].neighborSites[1] = boundary->neighborSites[(linkIdx + 1) % 3];
    vertices[currentCount].neighborSites[2] = boundary->neighborSites[(linkIdx + 2) % 3];
    
    /* Update connectivity */
    vertices[currentCount].adjacent[2] = inside;
    int k = 0;
    while (inside->adjacent[k] != boundary) k = (k + 1) % 3;
    inside->adjacent[k] = &vertices[currentCount];
    inside->neighborSites[k] = neighborIdx;
    
    Vertex3D *firstNew = &vertices[currentCount];
    Vertex3D *prevNew = firstNew;
    currentCount++;
    
    /* Continue around the cut boundary */
    Vertex3D *current = boundary;
    int currentLink = (linkIdx + 2) % 3;
    
    while (1) {
        Vertex3D *next = current->adjacent[currentLink];
        real_t offset;
        
        if (classifyVertexPosition3D(&neighbors[neighborIdx], next, weightFrac, &offset) == POINT_OUTSIDE) {
            /* Continue to next edge */
            int j = 0;
            while (next->adjacent[j] != current) j = (j + 1) % 3;
            currentLink = (j + 2) % 3;
            current = next;
        } else {
            /* Found inside vertex - create intersection */
            vertices[currentCount] = findEdgePlaneIntersection3D(current, next,
                                                                 &neighbors[neighborIdx], weightFrac);
            vertices[currentCount].status = STATUS_ACTIVE;
            vertices[currentCount].index = currentCount;
            
            /* Update connectivity */
            k = 0;
            while (next->adjacent[k] != current) k = (k + 1) % 3;
            
            vertices[currentCount].adjacent[2] = next;
            vertices[currentCount].neighborSites[2] = next->neighborSites[(k + 2) % 3];
            
            next->adjacent[k] = &vertices[currentCount];
            next->neighborSites[k] = neighborIdx;
            
            /* Link to previous new vertex */
            vertices[currentCount].adjacent[1] = prevNew;
            prevNew->adjacent[0] = &vertices[currentCount];
            
            vertices[currentCount].neighborSites[0] = neighborIdx;
            vertices[currentCount].neighborSites[1] = current->neighborSites[(currentLink + 1) % 3];
            
            prevNew = &vertices[currentCount];
            currentCount++;
            
            /* Continue traversal */
            int j = 0;
            while (next->adjacent[j] != &vertices[currentCount - 1]) j = (j + 1) % 3;
            currentLink = (j + 2) % 3;
            current = next;
            
            if (current == inside) break;
        }
        
        if (currentCount > 1000) {
            fprintf(stderr, "Warning: Vertex count exceeded limit\n");
            break;
        }
    }
    
    /* Close the loop */
    firstNew->adjacent[1] = prevNew;
    prevNew->adjacent[0] = firstNew;
    
    return currentCount;
}

/*============================================================================
 * MAIN CELL CONSTRUCTION - 3D
 *============================================================================*/

/**
 * @brief Construct 3D Laguerre cell for a given site
 */
int constructLaguerreCell3D(Point3D *center, Point3D *inputNeighbors,
                            Point3D *neighbors, int numNeighbors,
                            Vertex3D *vertices, int maxVertices,
                            real_t boxSize, int shrinkFlag, int *workArray) {
    (void)shrinkFlag;  /* Currently unused - reserved for future shrink operation */
    
    real_t maxDistSq = (boxSize/2) * (boxSize/2) * 3;
    real_t centerWeight = center->weightSquared;
    int vertexCount;
    
    /* Transform neighbors to local coordinates */
    int validCount = 0;
    real_t minDistSq = 1.0e20;
    
    for (int i = 0; i < numNeighbors; i++) {
        if (inputNeighbors[i].index != center->index) {
            neighbors[validCount] = subtract3D(&inputNeighbors[i], center);
            neighbors[validCount].soundSpeed = inputNeighbors[i].soundSpeed;
            neighbors[validCount].pressure = inputNeighbors[i].pressure;
            neighbors[validCount].squaredDist = dot3D(&neighbors[validCount], &neighbors[validCount]);
            minDistSq = MIN(minDistSq, neighbors[validCount].squaredDist);
            neighbors[validCount].weightSquared = inputNeighbors[i].weightSquared;
            
            real_t a = 0.5 * (1 + (centerWeight - neighbors[validCount].weightSquared) / neighbors[validCount].squaredDist);
            neighbors[validCount].deltaRadiusSq = neighbors[validCount].squaredDist * a * a;
            validCount++;
        }
    }
    numNeighbors = validCount;
    
    /* Sort by delta radius */
    qsort(neighbors, numNeighbors, sizeof(Point3D), compareByDeltaRadius3D);
    
    /* Initialize bounding box */
    vertexCount = initBoundingBox3D(vertices, maxVertices, boxSize);
    
    /* Cut cell with each neighbor */
    for (int i = 0; i < numNeighbors; i++) {
        if (neighbors[i].deltaRadiusSq > maxDistSq) break;
        
        real_t distSq = dot3D(&neighbors[i], &neighbors[i]);
        real_t neighborWeight = neighbors[i].weightSquared;
        real_t weightFrac = computeWeightFraction(centerWeight, neighborWeight, distSq);
        
        for (int j = 0; j < vertexCount; j++) {
            real_t offset;
            if (vertices[j].status == STATUS_ACTIVE &&
                classifyVertexPosition3D(&neighbors[i], &vertices[j], weightFrac, &offset) == POINT_OUTSIDE) {
                
                VertexRef3D start = findNearestBoundaryVertex(&neighbors[i], &vertices[j], weightFrac);
                inactivateOutsideVertices(&neighbors[i], vertices, vertexCount, weightFrac);
                vertexCount = createNewVertices3D(neighbors, i, &start, vertices, vertexCount, weightFrac);
                vertexCount = compactActiveVertices3D(vertices, vertexCount, &maxDistSq, workArray);
                break;
            }
        }
    }
    
    center->weightCeiling = computeMaxWeightCeiling3D(vertices, neighbors, vertexCount);
    return vertexCount;
}

/*============================================================================
 * UTILITY FUNCTIONS
 *============================================================================*/

real_t computeLaguerreWeight3D(real_t x, real_t y, real_t z, real_t smoothing) {
    real_t r = sqrt(x*x + y*y + z*z);
    real_t rh = r / smoothing;
    if (rh < 1.0) {
        /* Inside kernel */
    } else if (rh < 2.0) {
        /* Transition zone */
    }
    return 0;
}

int extractFaceVertices3D(Vertex3D *vertices, int numVertices, int faceId, Point3D *output) {
    int outputCount = 0;
    
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < 3; j++) {
            if (vertices[i].neighborSites[j] == faceId) {
                Vertex3D *start = &vertices[i];
                Vertex3D *current = start;
                int k = (j + 1) % 3;
                Vertex3D *target = vertices[i].adjacent[k];
                Vertex3D *startTarget = target;
                
                copyPoint3D(&output[outputCount], current);
                outputCount++;
                
                do {
                    for (k = 0; k < 3; k++) {
                        if (target->neighborSites[k] == faceId) {
                            current = target;
                            target = target->adjacent[(k + 1) % 3];
                            break;
                        }
                    }
                    copyPoint3D(&output[outputCount], current);
                    outputCount++;
                } while (current != start || target != startTarget);
                
                return outputCount;
            }
        }
    }
    
    return outputCount;
}

int extractPolygonVertices2D(Corner2D *start, int numCorners, Point2D *output) {
    Corner2D *current = start;
    
    for (int i = 0; i < numCorners; i++) {
        output[i].x = current->x;
        output[i].y = current->y;
        output[i].index = current->upperNeighborIdx;
        current = current->upperLink;
    }
    
    return numCorners;
}
