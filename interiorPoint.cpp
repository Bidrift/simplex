#include <bits/stdc++.h>

using namespace std;

class ColumnVector {
public:
    vector<double> columnValues;
    unsigned long long columnSize;
    ColumnVector() {
        columnValues = vector<double>();
        columnSize = 0;
    }
    ColumnVector(ColumnVector const &C) {
        columnValues = C.columnValues;
        columnSize = C.columnSize;
    }
    ColumnVector(ColumnVector &C) {
        columnSize = C.columnSize;
        columnValues = C.columnValues;
    }
    ColumnVector(vector<double>& columnValues) {
        this->columnValues = columnValues;
        columnSize = columnValues.size();
    }
    ColumnVector(unsigned long long size) {
        columnSize = size;
        columnValues = vector<double>(size);
    }
    ColumnVector operator+(ColumnVector& b) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] + b.columnValues[i];
        }
        return result;
    }
    ColumnVector operator-(ColumnVector& b) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] - b.columnValues[i];
        }
        return result;
    }
    ColumnVector operator*(double scalar) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] * scalar;
        }
        return result;
    }
    ColumnVector operator/(double scalar) {
        ColumnVector result(columnSize);
        for (int i = 0; i < columnSize; i++) {
            result.columnValues[i] = columnValues[i] / scalar;
        }
        return result;
    }
    friend istream& operator>>(istream& stream, ColumnVector& columnVector)
    {
        for (int i = 0; i < columnVector.columnSize; i++) {
            stream >> columnVector.columnValues[i];
        }
        return stream;
    }

    friend ostream& operator<<(ostream& stream, ColumnVector& columnVector)
    {
        for (int i = 0; i < columnVector.columnSize; i++) {
            stream << columnVector.columnValues[i] << endl;
        }
        return stream;
    }
    double norm() {
        double result = 0;
        for (int i = 0; i < columnSize; i++) {
            result += (columnValues[i]*columnValues[i]);
        }
        return sqrt(result);
    }
};


class Matrix {
protected:
    int nbRows, nbColumns;
    vector<ColumnVector> elements;
public:
    int getNbRows() const {
        return nbRows;
    }
    int getNbColumns() const {
        return nbColumns;
    }
    void setNbRows(int newNbRows) {
        nbRows = newNbRows;
    }
    void setNbColumns(int newNbColumns) {
        nbColumns = newNbColumns;
    }
    Matrix(int nbRows, int nbColumns) {
        this->nbRows = nbRows;
        this->nbColumns = nbColumns;
        for (int i = 0; i < nbColumns; i++) {
            elements.emplace_back(nbRows);
            for (int j = 0; j < nbColumns; j++){
                elements[i].columnValues.push_back(0);
            }
        }
    }
    Matrix(const ColumnVector& c){
        this->nbRows = c.columnSize;
        this->nbColumns = c.columnSize;
        for (int i = 0; i < nbColumns; i++) {
            elements.emplace_back(nbRows);
            for (int j = 0; j < nbColumns; j++){
                elements[i].columnValues.push_back(0);
            }
        }
        for (int i = 0; i < c.columnSize; i++) {
            this->setElement(i, i, c.columnValues[i]);
        }
    }
    void setElement(int i, int j, double value) {
        this->elements[j].columnValues[i] = value;
    }
    double& getElement(int i, int j) {
        return this->elements[j].columnValues[i];
    }
    Matrix operator+(Matrix& m1) {

        if (getNbRows() != m1.getNbRows() || getNbColumns() != m1.getNbColumns()) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix Empty(0,0);
            return Empty;
        }
        Matrix D(nbRows,nbColumns);
        for (int i = 0; i < nbRows; i++) {
            for (int j = 0; j < nbColumns; j++) {
                D.setElement(i,j, getElement(i,j) + m1.getElement(i,(j)));
            }
        }
        return D;
    }
    Matrix operator-(Matrix& m1) {
        if (getNbRows() != m1.getNbRows() || getNbColumns() != m1.getNbColumns()) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix Empty(0,0);
            return Empty;
        }
        Matrix D(nbRows,nbColumns);
        for (int i = 0; i < nbRows; i++) {
            for (int j = 0; j < nbColumns; j++) {
                D.setElement(i,j, getElement(i,j) - m1.getElement(i,(j)));
            }
        }
        return D;
    }
    Matrix operator*(Matrix& m1) {
        if (this->getNbColumns() != m1.getNbRows()) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix Empty(0,0);
            return Empty;
        }
        Matrix F(this->getNbRows(),m1.getNbColumns());
        for (int i = 0; i < this->getNbRows(); i++) {
            for (int j = 0; j < m1.getNbColumns(); j++) {
                double result = 0;
                for (int k = 0; k < nbColumns; k++){
                    result += this->getElement(i, k)*m1.getElement(k,j);
                }
                F.setElement(i,j,result);
            }
        }
        return F;
    }

    ColumnVector operator*(ColumnVector& C) {
        if (this->getNbColumns() != C.columnSize) {
            cout << "Error: the dimensional problem occurred\n";
            ColumnVector Empty(0);
            return Empty;
        }
        ColumnVector F(getNbRows());
        for (int i = 0; i < getNbRows(); i++) {
            double result = 0;
            for (int j = 0; j < C.columnSize; j++) {
                result += getElement(i,j)*C.columnValues[j];
            }
            F.columnValues[i] = result;
        }
        return F;
    }

    Matrix& operator=(Matrix m1) {
        this->setNbRows(m1.getNbRows());
        this->setNbColumns(m1.getNbColumns());
        for (int i = 0; i < nbRows; i++) for (int j = 0; j < nbColumns; j++) {
                this->setElement(i,j, m1.getElement(i,j));
            }
        return *this;
    }

    friend istream& operator>>(istream& stream, Matrix& matrix)
    {
        for (int i = 0; i < matrix.getNbRows(); i++) for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream >> matrix.getElement(i,j);
            }
        return stream;
    }

    friend ostream& operator<<(ostream& stream, Matrix& matrix)
    {
        for (int i = 0; i < matrix.getNbRows(); i++) {
            for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream << matrix.getElement(i,j) << " ";
            }
            stream << endl;
        }
        return stream;
    }
    Matrix transpose() {
        Matrix G(nbColumns,nbRows);
        for (int i = 0; i < nbColumns; i++) {
            for (int j = 0; j < nbRows; j++) {
                G.setElement(i,j, getElement(j,i));
            }
        }
        return G;
    }
};

class SquareMatrix: public Matrix {
public:
    explicit SquareMatrix(int size) : Matrix(size, size) {};
    SquareMatrix(Matrix matrix) : Matrix(matrix.getNbRows(), matrix.getNbColumns()) {
        for (int i = 0; i < nbRows; i++)
            for (int j = 0; j < nbColumns; j++)
                setElement(i, j, matrix.getElement(i, j));
    }
    SquareMatrix operator+(SquareMatrix& sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        Matrix result = *m + *m1;
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
    SquareMatrix operator-(SquareMatrix& sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        Matrix result = *m - *m1;
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
    SquareMatrix operator*(SquareMatrix& sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        Matrix result = *m * *m1;
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
    ColumnVector operator*(ColumnVector& C) {
        auto *m = (Matrix*) this;
        ColumnVector *c = &C;
        ColumnVector result = (Matrix)(*m) * *c;
        auto *sResult = (ColumnVector *) &result;
        return *sResult;
    }
    SquareMatrix& operator=(SquareMatrix sm1) {
        Matrix *m = this;
        Matrix *m1 = &sm1;
        *m = *m1;
        return *this;
    }
    SquareMatrix transpose() {
        Matrix *m = this;
        Matrix result = m->transpose();
        auto* sResult = (SquareMatrix *) &result;
        return *sResult;
    }
};

class IdentityMatrix: public SquareMatrix {
public:
    explicit IdentityMatrix(int size): SquareMatrix(size) {
        for (int i = 0; i < getNbRows(); i++) for (int j = 0; j < getNbColumns(); j++) {
                if (i == j) {
                    this->setElement(i,j,1);
                } else {
                    this->setElement(i,j,0);
                }
            }
    };
};

class EliminationMatrix: public SquareMatrix {
private:

    int iNullified;
    int jNullified;
public:
    EliminationMatrix(SquareMatrix &matrix, int i, int j) : SquareMatrix(matrix.getNbRows()) {
        IdentityMatrix identityMatrix(matrix.getNbRows());
        *this = *((EliminationMatrix *)(SquareMatrix *) (&identityMatrix));
        iNullified = i;
        jNullified = j;
        this->setElement(i,j,-(matrix.getElement(i,j)/matrix.getElement(j,j)));
    }
    int getINullified() const {
        return iNullified;
    }

    int getJNullified() const {
        return jNullified;
    }

};

class PermutationMatrix: public SquareMatrix {
private:
    int firstSwapped;
    int secondSwapped;
    void swapRows(int i1, int i2) {
        for (int i = 0; i < nbRows; i++) {
            swap(getElement(i1,i), getElement(i2,i));
        }
    }
public:
    PermutationMatrix(SquareMatrix& matrix, int i1, int i2) : SquareMatrix(matrix.getNbRows()) {
        IdentityMatrix identityMatrix(matrix.getNbRows());
        *this = *((PermutationMatrix *)(SquareMatrix *) (&identityMatrix));
        firstSwapped = i1;
        secondSwapped = i2;
        swapRows(firstSwapped,secondSwapped);
    }
};

class AugmentedMatrix: public Matrix {
public:
    SquareMatrix m1;
    SquareMatrix m2;
    explicit AugmentedMatrix(SquareMatrix &M1) : Matrix(M1.getNbRows(), M1.getNbColumns()), m1(M1), m2(IdentityMatrix(M1.getNbRows())) {};

    friend ostream& operator<<(ostream& stream, AugmentedMatrix& matrix)
    {
        for (int i = 0; i < matrix.getNbRows(); i++) {
            for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream << fixed << setprecision(2) << matrix.m1.getElement(i,j) << " ";
            }
            for (int j = 0; j < matrix.getNbColumns(); j++) {
                stream << fixed << setprecision(2) << matrix.m2.getElement(i,j) << " ";
            }
            stream << endl;
        }
        return stream;
    }
};

class FindInverseMatrix {
private:
    int step;

    static int findPermutation(AugmentedMatrix& matrix, int column) {
        int pivot = column;
        for (int j = column+1; j < matrix.getNbColumns(); j++) {
            if (abs(matrix.m1.getElement(j,column)) > abs(matrix.m1.getElement(pivot,column))) {
                pivot = j;
            }
        }
        return pivot;
    }
    void makePermutation(AugmentedMatrix& matrix, int row1, int row2) {
        if (row1 != row2) {
            PermutationMatrix P(matrix.m1, row1, row2);
            matrix.m1 = P*matrix.m1;
            matrix.m2 = P*matrix.m2;
            step++;
        }
    }
    void eliminatePosition(AugmentedMatrix& matrix, int i, int j) {
        EliminationMatrix E(matrix.m1,i,j);
        matrix.m1 = E*matrix.m1;
        matrix.m2 = E*matrix.m2;
        step++;
    }

    void directWay(AugmentedMatrix& augmentedMatrix) {
        for (int i = 0; i < augmentedMatrix.getNbColumns()-1; i++) {
            int pivot = findPermutation(augmentedMatrix,i);
            makePermutation(augmentedMatrix,i,pivot);
            for (int j = i+1; j < augmentedMatrix.getNbColumns(); j++) {
                if (augmentedMatrix.m1.getElement(j,i) != 0) {
                    eliminatePosition(augmentedMatrix,j,i);
                }
            }
        }
    }

    void wayBack(AugmentedMatrix& augmentedMatrix) {
        for (int i = augmentedMatrix.getNbColumns()-1; i > 0; i--) {
            for (int j = i-1; j >= 0; j--) {
                if (augmentedMatrix.m1.getElement(j,i) != 0) {
                    eliminatePosition(augmentedMatrix,j,i);
                }
            }
        }
    }

    static void normalizeDiagonal(AugmentedMatrix& augmentedMatrix) {
        for (int i = 0; i < augmentedMatrix.getNbRows(); i++) {
            double pivot = augmentedMatrix.m1.getElement(i,i);
            augmentedMatrix.m1.setElement(i,i,augmentedMatrix.m1.getElement(i,i)/pivot);
            for (int j = 0; j < augmentedMatrix.getNbColumns(); j++) {
                augmentedMatrix.m2.setElement(i,j,augmentedMatrix.m2.getElement(i,j)/pivot);
            }
        }
    }
public:
    SquareMatrix findInverse(SquareMatrix matrix) {
        step = 0;
        AugmentedMatrix augmentedMatrix(matrix);
        step++;
        directWay(augmentedMatrix);
        wayBack(augmentedMatrix);
        normalizeDiagonal(augmentedMatrix);
        matrix = augmentedMatrix.m2;
        return matrix;
    }
};


void solve(ColumnVector c, Matrix A, ColumnVector rightHand, ColumnVector initial, double alpha, int nb_variables, int nb_equations, int nb_slack, int itr) {
    Matrix D(initial);
    Matrix A_ = A*D;
    ColumnVector c_ = D*c;
    Matrix A_T = A_.transpose();
    SquareMatrix A_A_T(A_*A_T);
    FindInverseMatrix inverseMatrix;
    Matrix A_A_T_1 = inverseMatrix.findInverse(A_A_T);
    IdentityMatrix I(nb_variables+nb_slack);
    Matrix P = A_A_T_1*A_;
    P = A_T*P;
    P = Matrix(I) - P;
    ColumnVector c_p = P*c_;
    double mini = 0;
    for (int i = 0; i < nb_slack + nb_variables; i++) {
        mini = min(c_p.columnValues[i],mini);
    }
    double v = abs(mini);
    ColumnVector ones(nb_variables + nb_slack);
    for (int i = 0; i < nb_variables + nb_slack; i++) {
        ones.columnValues[i] = 1;
    }
    ColumnVector x_ = c_p*(alpha/v);
    x_ = x_ + ones;
    ColumnVector x = D*x_;
    cout << "Iteration #" << itr << endl;
    cout << x;
    bool found = true;
    for (int i = 0; i < nb_slack + nb_variables; i++) {
        if (abs(x.columnValues[i] - initial.columnValues[i]) > 1e-4)
            found = false;
    }
    if (found) {
        return;
    } else {
        solve(c, A, rightHand, x, alpha, nb_variables, nb_equations, nb_slack, itr + 1);
    }
}

bool checkValid(Matrix A, ColumnVector rightHand, ColumnVector initial, int nb_variables, int nb_equations, int nb_slack) {
    for (int i = 0; i < nb_variables+nb_slack; i++) {
        if (initial.columnValues[i] <= 0) return false;
    }
    for (int i = 0; i < nb_equations; i++) {
        double sum = 0;
        for (int j = 0; j < nb_variables + nb_slack; j++) {
            sum += A.getElement(i, j)*initial.columnValues[j];
        }
        if (sum != rightHand.columnValues[i]) return false;
    }
    return true;
}


int main() {
    int minimize = 2;
    while (minimize != 1 && minimize != 0) {
        cout << "Please specify whether this is a minimization or maximization problem.\n1 - Minimization\n0 - Maximization\n";
        cin >> minimize;
    }
    int nb_variables = 0;
    while (nb_variables <= 0) {
        cout << "Please specify the number of variables\n";
        cin >> nb_variables;
    }
    cout << "Please submit the coefficients for the z expression\n";
    ColumnVector coefficient_c(nb_variables);
    cin >> coefficient_c;
    if (minimize) {
        coefficient_c = coefficient_c*(-1);
    }
    int nb_equations = 0;
    while (nb_equations <= 0) {
        cout << "Please specify the number of equations (don't count sign constraints)\n";
        cin >> nb_equations;
    }
    vector<int> slack;
    int nb_slack = 0;
    for (int i = 0; i < nb_equations; i++) {
        cout << "Please specify the type of inequation for inequation number " << i+1 << "\n0 - equal\n1 - less\n2 - greater\n";
        int x; cin >> x;
        while (x < 0 || x > 2) {
            cout << "Please specify the type of inequation for inequation number " << i+1 << "\n0 - equal\n1 - less\n2 - greater\n";
            cin >> x;
        }
        if (x == 0) slack.push_back(0);
        else if (x == 1) {
            nb_slack++;
            slack.push_back(1);
        } else {
            nb_slack++;
            slack.push_back(-1);
        }
    }
    Matrix A_low(nb_equations, nb_variables);
    for (int i = 0; i < nb_equations; i++) {
        cout << "Equation " << i+1 << "\n";
        for (int j = 0; j < nb_variables; j++) {
            int x; cin >> x;
            A_low.setElement(i, j, x);
        }
    }
    ColumnVector rightHand(nb_equations);
    cout << "Please specify the right hand values\n";
    cin >> rightHand;
    double alpha;
    cout << "Please specify the value for alpha\n";
    cin >> alpha;
    Matrix A(nb_equations, nb_variables + nb_slack);
    for (int i = 0; i < nb_equations; i++)
        for (int j = 0; j < nb_variables; j++)
            A.setElement(i, j, A_low.getElement(i, j));
    int j = 0;
    for (int i = 0; i < slack.size(); i++) {
        A.setElement(i, j+nb_variables, slack[i]);
        if (slack[i] != 0) j++;
    }
    ColumnVector c(nb_variables+nb_slack);
    for (int i = 0; i < nb_variables; i++) {
        c.columnValues[i] = coefficient_c.columnValues[i];
    }
    ColumnVector initial(nb_variables + nb_slack);
    while (!checkValid(A, rightHand, initial, nb_variables, nb_equations, nb_slack)) {
        cout << "Please specify the initial values for the Interior-Point algorithm\n";
        cin >> initial;
    }
    solve(c, A, rightHand, initial, alpha, nb_variables, nb_equations, nb_slack, 1);
}
