/* MATRIX INDEXING, BY FAR OUR BEST OPTION FOR BIGGER MATRIXES */
/* EXECUTION TIMES
    test    time
    0       0.010s
    1       0.378s
    2       0.415s
    30      8.187s
    500     ~116m
    1000    19m47.414s
*/

#include <stdlib.h>
#include <bits/stdc++.h>

#define RAND01 ((double)random() / (double)RAND_MAX)
#define TRUE 1
#define FALSE 0

using namespace std;

typedef double** matrix;
typedef matrix* matrixPtr;
typedef vector<int> vi; 
typedef vector<double> vd; 
typedef struct {
    int i;
    int j;
    int val;
} s_entry;

typedef vector<s_entry> s_matrix;


matrix L, R;
vector<s_matrix> entries_in_row;
vector<s_matrix> entries_in_column;

void initialize_matrix(matrixPtr M, int lines, int columns);
void random_fill_LR(int nU, int nI, int nF);
void copy(matrix M, matrix M_aux, int lines, int columns);
void transpose(matrix M, matrix M_t, int lines, int columns);
void product_transposed(matrix R_t, matrix B, int nU, int nI, int nF);
void iterate(int numIter, double alpha, int nU, int nI, int nF, matrix L_aux, matrix R_aux, matrix R_t, s_matrix A, matrix B);
void print_matrix(matrix M, int lines, int columns);
void print_s_matrix(s_matrix A);
void print_result(s_matrix A, matrix B, int lines, int columns);
void free_matrix(matrix M, int lines);
int coord_in_A(s_matrix A, int line, int column);

// Allocate memory for MxN matrixes
void initialize_matrix(matrixPtr M, int lines, int columns)
{
    *M = (double **) malloc(lines * sizeof(double *)); 
    for (int i = 0; i < lines; i++) 
        (*M)[i] = (double *) malloc(columns * sizeof(double)); 
}

// As provided in the statement
void random_fill_LR(int nU, int nI, int nF)
{
    srandom(0);

    for(int i = 0; i < nU; i++)
        for(int j = 0; j < nF; j++)
            L[i][j] = RAND01 / (double) nF;
    
    for(int i = 0; i < nF; i++) 
        for(int j = 0; j < nI; j++)
            R[i][j] = RAND01 / (double) nF;

}

void copy(matrix M, matrix M_aux, int lines, int columns)
{
    for (int i = 0; i < lines; i++)
        for (int j = 0; j < columns; j++)
            M_aux[i][j] = M[i][j];
}

void transpose(matrix M, matrix M_t, int lines, int columns)
{
    for (int i = 0; i < lines; i++)
        for (int j = 0; j < columns; j++)
            M_t[j][i] = M[i][j];
}

// Multiply L x R_t
// Using R_t instead of R enhances hit rates in cache
void product_transposed(matrix R_t, matrix B, int nU, int nI, int nF)
{
    double sum = 0;
    for (int i = 0; i < nU; i++) {
        for (int j = 0; j < nI; j++) {
            for (int k = 0; k < nF; k++) {
                sum += L[i][k] * R_t[j][k];
            }
            B[i][j] = sum;
            sum = 0;
        }
    }
}

void print_s_matrix(s_matrix A) 
{ 
    for_each(A.begin(), A.end(), [](s_entry a) {
        cout << a.i << " " << a.j << " " << a.val << endl;
    });
} 

void print_matrix(matrix M, int lines, int columns) 
{ 
    for (int i = 0; i < lines; i++) { 
        for (int j = 0; j < columns; j++)  
            cout << M[i][j] << " ";         
        cout << endl; 
    } 
}

// check if value is in vector A
int coord_in_A(s_matrix A, int line, int column) {
    for (int p = 0; p < A.size(); p++) {
        if (A[p].i == line && A[p].j == column)
            return TRUE;
        else if ((A[p].i == line && A[p].j > column) || (A[p].i > line)) {
            return FALSE;
        }
    }
    return FALSE;
}

// index matrix A by row: store all values in a single row
s_matrix indexes_of_row(s_matrix A, int line) {
    s_matrix result;
    for (int p = 0; p < A.size(); p++) {
        if (A[p].i == line) 
            result.push_back(A[p]);
        else if (A[p].i > line) // assumption that input file has matrix items ordered by line
            break;
    }
    return result;
}

// index matrix A by column: store all values in a single column
s_matrix indexes_of_column(s_matrix A, int column) {
    s_matrix result;
    for (int p = 0; p < A.size(); p++)
        if (A[p].j == column) 
            result.push_back(A[p]);
    
    return result;
}

// index matrix to simplify inner for cycle in iterate
void index_matrix(s_matrix A, int lines, int columns)
{
    entries_in_row.resize(lines);
    entries_in_column.resize(columns);
    for (int p = 0; p < A.size(); p++) {

        if (entries_in_row[A[p].i].empty())
            entries_in_row[A[p].i] = indexes_of_row(A, A[p].i);

        if (entries_in_column[A[p].j].empty())
            entries_in_column[A[p].j] = indexes_of_column(A, A[p].j);

    }
}

// free memory allocated
void free_matrix(matrix M, int lines)
{
    for (int i = 0; i < lines; i++)
        free(M[i]);
    free(M);
}

/* Strassen's Algorithm for matrix multiplication 
Complexity: O(n^2.808) */

inline double** MatrixMultiply(double** a, double** b, int n, 
									int l, int m) 
{ 
	double** c = new double*[n]; 
	for (int i = 0; i < n; i++) 
		c[i] = new double[m]; 

	for (int i = 0; i < n; i++) { 
		for (int j = 0; j < m; j++) { 
			c[i][j] = 0; 
			for (int k = 0; k < l; k++) { 
				c[i][j] += a[i][k] * b[k][j]; 
			} 
		} 
	} 
	return c; 
} 

inline double** Strassen(double** a, double** b, int n, 
								int l, int m) 
{ 
	if (n == 1 || l == 1 || m == 1) 
		return MatrixMultiply(a, b, n, l, m); 

	double** c = new double*[n]; 
	for (int i = 0; i < n; i++) 
		c[i] = new double[m]; 

	int adjN = (n >> 1) + (n & 1); 
	int adjL = (l >> 1) + (l & 1); 
	int adjM = (m >> 1) + (m & 1); 


	double**** As = new double***[2]; 
	for (int x = 0; x < 2; x++) { 
		As[x] = new double**[2]; 
		for (int y = 0; y < 2; y++) { 
			As[x][y] = new double*[adjN]; 
			for (int i = 0; i < adjN; i++) { 
				As[x][y][i] = new double[adjL]; 
				for (int j = 0; j < adjL; j++) { 
					int I = i + (x & 1) * adjN; 
					int J = j + (y & 1) * adjL; 
					As[x][y][i][j] = (I < n && J < l) ? a[I][J] : 0; 
				} 
			} 
		} 
	} 

	double**** Bs = new double***[2]; 
	for (int x = 0; x < 2; x++) { 
		Bs[x] = new double**[2]; 
		for (int y = 0; y < 2; y++) { 
			Bs[x][y] = new double*[adjN]; 
			for (int i = 0; i < adjL; i++) { 
				Bs[x][y][i] = new double[adjM]; 
				for (int j = 0; j < adjM; j++) { 
					int I = i + (x & 1) * adjL; 
					int J = j + (y & 1) * adjM; 
					Bs[x][y][i][j] = (I < l && J < m) ? b[I][J] : 0; 
				} 
			} 
		} 
	} 

	double*** s = new double**[10]; 
	for (int i = 0; i < 10; i++) { 
		switch (i) { 
		case 0: 
			s[i] = new double*[adjL]; 
			for (int j = 0; j < adjL; j++) { 
				s[i][j] = new double[adjM]; 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[0][1][j][k] - Bs[1][1][j][k]; 
				} 
			} 
			break; 
		case 1: 
			s[i] = new double*[adjN]; 
			for (int j = 0; j < adjN; j++) { 
				s[i][j] = new double[adjL]; 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][0][j][k] + As[0][1][j][k]; 
				} 
			} 
			break; 
		case 2: 
			s[i] = new double*[adjN]; 
			for (int j = 0; j < adjN; j++) { 
				s[i][j] = new double[adjL]; 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[1][0][j][k] + As[1][1][j][k]; 
				} 
			} 
			break; 
		case 3: 
			s[i] = new double*[adjL]; 
			for (int j = 0; j < adjL; j++) { 
				s[i][j] = new double[adjM]; 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[1][0][j][k] - Bs[0][0][j][k]; 
				} 
			} 
			break; 
		case 4: 
			s[i] = new double*[adjN]; 
			for (int j = 0; j < adjN; j++) { 
				s[i][j] = new double[adjL]; 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][0][j][k] + As[1][1][j][k]; 
				} 
			} 
			break; 
		case 5: 
			s[i] = new double*[adjL]; 
			for (int j = 0; j < adjL; j++) { 
				s[i][j] = new double[adjM]; 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[0][0][j][k] + Bs[1][1][j][k]; 
				} 
			} 
			break; 
		case 6: 
			s[i] = new double*[adjN]; 
			for (int j = 0; j < adjN; j++) { 
				s[i][j] = new double[adjL]; 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][1][j][k] - As[1][1][j][k]; 
				} 
			} 
			break; 
		case 7: 
			s[i] = new double*[adjL]; 
			for (int j = 0; j < adjL; j++) { 
				s[i][j] = new double[adjM]; 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[1][0][j][k] + Bs[1][1][j][k]; 
				} 
			} 
			break; 
		case 8: 
			s[i] = new double*[adjN]; 
			for (int j = 0; j < adjN; j++) { 
				s[i][j] = new double[adjL]; 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][0][j][k] - As[1][0][j][k]; 
				} 
			} 
			break; 
		case 9: 
			s[i] = new double*[adjL]; 
			for (int j = 0; j < adjL; j++) { 
				s[i][j] = new double[adjM]; 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[0][0][j][k] + Bs[0][1][j][k]; 
				} 
			} 
			break; 
		} 
	} 

	double*** p = new double**[7];
 
	p[0] = Strassen(As[0][0], s[0], adjN, adjL, adjM); 
	p[1] = Strassen(s[1], Bs[1][1], adjN, adjL, adjM); 
	p[2] = Strassen(s[2], Bs[0][0], adjN, adjL, adjM); 
	p[3] = Strassen(As[1][1], s[3], adjN, adjL, adjM); 
	p[4] = Strassen(s[4], s[5], adjN, adjL, adjM); 
	p[5] = Strassen(s[6], s[7], adjN, adjL, adjM); 
	p[6] = Strassen(s[8], s[9], adjN, adjL, adjM); 

	for (int i = 0; i < adjN; i++) { 
		for (int j = 0; j < adjM; j++) { 
			c[i][j] = p[4][i][j] + p[3][i][j] - p[1][i][j] + p[5][i][j]; 
			if (j + adjM < m) 
				c[i][j + adjM] = p[0][i][j] + p[1][i][j]; 
			if (i + adjN < n) 
				c[i + adjN][j] = p[2][i][j] + p[3][i][j]; 
			if (i + adjN < n && j + adjM < m) 
				c[i + adjN][j + adjM] = p[4][i][j] + p[0][i][j] - p[2][i][j] - p[6][i][j]; 
		} 
	} 

	for (int x = 0; x < 2; x++) { 
		for (int y = 0; y < 2; y++) { 
			for (int i = 0; i < adjN; i++) { 
				delete[] As[x][y][i]; 
			} 
			delete[] As[x][y]; 
		} 
		delete[] As[x]; 
	} 
	delete[] As; 

	for (int x = 0; x < 2; x++) { 
		for (int y = 0; y < 2; y++) { 
			for (int i = 0; i < adjL; i++) { 
				delete[] Bs[x][y][i]; 
			} 
			delete[] Bs[x][y]; 
		} 
		delete[] Bs[x]; 
	} 
	delete[] Bs; 

	for (int i = 0; i < 10; i++) { 
		switch (i) { 
		case 0: 
		case 3: 
		case 5: 
		case 7: 
		case 9: 
			for (int j = 0; j < adjL; j++) { 
				delete[] s[i][j]; 
			} 
			break; 
		case 1: 
		case 2: 
		case 4: 
		case 6: 
		case 8: 
			for (int j = 0; j < adjN; j++) { 
				delete[] s[i][j]; 
			} 
			break; 
		} 
		delete[] s[i]; 
	} 
	delete[] s; 

	for (int i = 0; i < 7; i++) { 
		for (int j = 0; j < (n >> 1); j++) { 
			delete[] p[i][j]; 
		} 
		delete[] p[i]; 
	} 
	delete[] p; 

	return c; 
} 

void iterate(int numIter, double alpha, int nU, int nI, int nF, matrix L_aux, matrix R_aux, matrix R_t, s_matrix A, matrix B)
{
    matrix tmp;
    // transpose(R, R_t, nF, nI);
    // product_transposed(R_t, B, nU, nI, nF);
    B=Strassen(L,R,nU,nF,nI);


    for (int count = 0; count < numIter; count++) {

        // auto start = chrono::high_resolution_clock::now();

        for (int p = 0; p < A.size(); p++) {
            int i = A[p].i, j = A[p].j;

            for (int k = 0; k < nF; k++) {
                double deltadl = 0, deltadr = 0;

                for (s_entry Aij : entries_in_row[i]) {
                    deltadl += 2 * (Aij.val - B[i][Aij.j]) * (-1) * (R[k][Aij.j]);
                }
                L_aux[i][k] = L[i][k] - alpha * deltadl;

                for (s_entry Aij : entries_in_column[j]) {
                    deltadr += 2 * (Aij.val - B[Aij.i][j]) * (-1) * (L[Aij.i][k]);              
                }
                R_aux[k][j] = R[k][j] - alpha * deltadr;
            } 
        }
        
        tmp = L;
        L = L_aux;
        L_aux = tmp;

        tmp = R;
        R = R_aux;
        R_aux = tmp;

        // transpose(R, R_t, nF, nI);
        // product_transposed(R_t, B, nU, nI, nF); 

        B=Strassen(L,R,nU,nF,nI);
  

        // auto finish = chrono::high_resolution_clock::now();
        // chrono::duration<double> elapsed = finish - start;
        // cout << "Finished iteration " << count << " in " << elapsed.count() << endl;
    }
}

void print_result(s_matrix A, matrix B, int lines, int columns)
{
    print_matrix(B, lines, columns);
    for(int i = 0; i < lines; i++) {
        double max = 0;
        int max_j = 0;
        for(int j = 0; j < columns; j++) {
            if (B[i][j] > max) {
                if (coord_in_A(A, i, j) == FALSE) {
                    max = B[i][j];
                    max_j = j;
                }
            }
        }
        cout << max_j << endl;
    }
}



int main(int argc, char** argv)
{

    string inFile = "";
    string outFile = "";

    int numIter = 0, nF = 0, nU = 0, nI = 0, NNZ = 0;
    int nr = 0, nc = 0;
    double alpha = 0.000, val = 0.000;

    s_matrix A;
    
    ifstream infile;

    if( argc == 2 ) {
      inFile = argv[1];
    }
    else {
      cout << "Usage: ./cfile InputFile\n";
      return 1;
    }

    infile.open(inFile);
    infile >> numIter;
    infile >> alpha;
    infile >> nF;
    infile >> nU >> nI >> NNZ;

    while(infile >> nr >> nc >> val) {
        s_entry tmp;
        tmp.i = nr;
        tmp.j = nc;
        tmp.val = val;
        A.push_back(tmp);
    }     

    infile.close();
    
    /* INITIALIZATION */

    auto start = chrono::high_resolution_clock::now();

    initialize_matrix(&L, nU, nF);
    initialize_matrix(&R, nF, nI);
    random_fill_LR(nU, nI, nF);

    index_matrix(A, nU, nI);


    matrix L_aux, R_aux, R_t, B;
    initialize_matrix(&L_aux, nU, nF);
    initialize_matrix(&R_aux, nF, nI);
    initialize_matrix(&R_t, nI, nF);
    initialize_matrix(&B, nU, nI);

    copy(L, L_aux, nU, nF);
    copy(R, R_aux, nF, nI);

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Seconds in preparation:" << elapsed.count() << endl;

    /* PROCESSING */

    start = chrono::high_resolution_clock::now();

    iterate(numIter, alpha, nU, nI, nF, L_aux, R_aux, R_t, A, B);

    finish = chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "Seconds in iteration:" << elapsed.count() << endl;

    /* RESULT */

    print_result(A, B, nU, nI);

    /* CLEANUP */

    free_matrix(L, nU);
    free_matrix(R, nF);
    free_matrix(L_aux, nU);
    free_matrix(R_aux, nF);
    free_matrix(R_t, nI);
    free_matrix(B, nU);

    return 0;
}