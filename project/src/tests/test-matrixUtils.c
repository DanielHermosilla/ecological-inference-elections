#include "../utils/matrixUtils.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void test_makeArray()
{
    double array[5];
    makeArray(array, 5, 3.14);
    for (int i = 0; i < 5; i++)
    {
        assert(array[i] == 3.14);
    }
    printf("makeArray() test passed.\n");
}

void test_createMatrix()
{
    Matrix m = createMatrix(3, 3);
    assert(m.rows == 3);
    assert(m.cols == 3);
    assert(m.data != NULL);
    freeMatrix(&m);
    printf("createMatrix() test passed.\n");
}

void test_checkMatrix()
{
    Matrix m = createMatrix(2, 2);
    checkMatrix(&m); // Should not throw an error
    freeMatrix(&m);
    printf("checkMatrix() test passed.\n");
}

void test_fillMatrix()
{
    Matrix m = createMatrix(2, 2);
    fillMatrix(&m, 5.0);
    for (int i = 0; i < m.rows * m.cols; i++)
    {
        assert(m.data[i] == 5.0);
    }
    freeMatrix(&m);
    printf("fillMatrix test passed.\n");
}

void test_rowSum()
{
    Matrix m = createMatrix(2, 3);
    fillMatrix(&m, 1.0);
    double result[2];
    rowSum(&m, result);
    assert(result[0] == 3.0);
    assert(result[1] == 3.0);
    freeMatrix(&m);
    printf("rowSum() test passed.\n");
}

void test_colSum()
{
    Matrix m = createMatrix(2, 3);
    fillMatrix(&m, 2.0);
    double result[3];
    colSum(&m, result);
    for (int i = 0; i < 3; i++)
    {
        assert(result[i] == 4.0);
    }
    freeMatrix(&m);
    printf("colSum() test passed.\n");
}

void test_maxElement()
{
    Matrix m = createMatrix(2, 2);
    m.data[0] = 1.0;
    m.data[1] = 5.0;
    m.data[2] = 3.0;
    m.data[3] = 2.0;
    double max = maxElement(&m);
    assert(max == 5.0);
    freeMatrix(&m);
    printf("maxElement() test passed.\n");
}

void test_convergeMatrix()
{
    Matrix m1 = createMatrix(2, 2);
    Matrix m2 = createMatrix(2, 2);
    fillMatrix(&m1, 1.0);
    fillMatrix(&m2, 1.01);
    bool converges = convergeMatrix(&m1, &m2, 0.02);
    bool converges2 = convergeMatrix(&m1, &m2, 0.001);
    bool converges3 = convergeMatrix(&m1, &m2, 0.01);

    assert(converges);
    assert(!converges2);
    assert(!converges3); // The converge criteria is an <=

    freeMatrix(&m1);
    freeMatrix(&m2);
    printf("convergeMatrix() test passed.\n");
}

void test_removeLastRow()
{
    Matrix m = createMatrix(3, 3);
    fillMatrix(&m, 1.0);
    Matrix reduced = removeLastRow(&m);
    assert(reduced.rows == 2);
    assert(reduced.cols == 3);
    freeMatrix(&m);
    freeMatrix(&reduced);
    printf("removeLastRow() test passed.\n");
}

void test_removeLastColumn()
{
    Matrix m = createMatrix(3, 3);
    fillMatrix(&m, 1.0);
    Matrix reduced = removeLastColumn(&m);
    assert(reduced.rows == 3);
    assert(reduced.cols == 2);
    freeMatrix(&m);
    freeMatrix(&reduced);
    printf("removeLastColumn() test passed.\n");
}

void test_createDiagonalMatrix()
{
    double values[3] = {1.0, 2.0, 3.0};
    Matrix diag = createDiagonalMatrix(values, 3);
    assert(diag.rows == 3);
    assert(diag.cols == 3);
    assert(diag.data[0] == 1.0 && diag.data[4] == 2.0 && diag.data[8] == 3.0);
    freeMatrix(&diag);
    printf("createDiagonalMatrix() test passed.\n");
}

void test_copyMatrix()
{
    Matrix m = createMatrix(2, 2);
    fillMatrix(&m, 3.0);
    Matrix copy = copyMatrix(&m);
    assert(copy.rows == m.rows);
    assert(copy.cols == m.cols);
    for (int i = 0; i < m.rows * m.cols; i++)
    {
        assert(copy.data[i] == m.data[i]);
    }
    freeMatrix(&m);
    freeMatrix(&copy);
    printf("copyMatrix() test passed.\n");
}

void test_getRow()
{
    Matrix m = createMatrix(2, 2);
    m.data[0] = 1.0;
    m.data[1] = 2.0;
    m.data[2] = 3.0;
    m.data[3] = 4.0;
    double *row = getRow(&m, 1);
    assert(row[0] == 3.0 && row[1] == 4.0);
    free(row);
    freeMatrix(&m);
    printf("getRow() test passed.\n");
}

void test_getColumn()
{
    Matrix m = createMatrix(2, 2);
    m.data[0] = 1.0;
    m.data[1] = 2.0;
    m.data[2] = 3.0;
    m.data[3] = 4.0;
    double *column = getColumn(&m, 1);
    assert(column[0] == 2.0 && column[1] == 4.0);
    free(column);
    freeMatrix(&m);
}

void test_addRowToMatrix()
{
    Matrix m = createMatrix(2, 3); // 2x3 matrix
    double newRow[3] = {4.0, 5.0, 6.0};
    addRowToMatrix(&m, newRow);
    assert(m.rows == 3);
    for (int j = 0; j < m.cols; j++)
    {
        assert(MATRIX_AT(m, 2, j) == newRow[j]);
    }
    freeMatrix(&m);
}

void test_addRowOfZeros()
{
    Matrix m = createMatrix(2, 2);
    fillMatrix(&m, 3.0);
    addRowOfZeros(&m, 1);
    for (int j = 0; j < m.cols; j++)
    {
        assert(MATRIX_AT(m, 1, j) == 0);
    }
    freeMatrix(&m);
}

void test_removeRow()
{
    Matrix m = createMatrix(8, 3); // 2x3 matrix
    for (int j = 0; j < m.cols; j++)
    {
        MATRIX_AT(m, 4, j) = 10;
    }
    removeRow(&m, 4);
    assert(m.rows == 7);
    for (int j = 0; j < m.cols; j++)
    {
        assert(MATRIX_AT(m, 4, j) != 10);
    }
    freeMatrix(&m);
}

void test_addColumnOfZeros()
{
    Matrix m = createMatrix(8, 3); // 2x3 matrix
    addColumnOfZeros(&m, 3);
    assert(m.cols == 4);
    freeMatrix(&m);
}

void test_removeColumn()
{
    Matrix m = createMatrix(8, 3); // 2x3 matrix
    for (int k = 0; k < m.rows; k++)
    {
        MATRIX_AT(m, k, 1) = 10;
    }
    removeColumn(&m, 1);
    assert(m.cols == 2);
    for (int k = 0; k < m.rows; k++)
    {
        MATRIX_AT(m, k, 1) = 0;
    }
    freeMatrix(&m);
}

int main()
{
    printf("Running matrixUtils tests...\n");
    test_makeArray();
    test_createMatrix();
    test_checkMatrix();
    test_fillMatrix();
    test_rowSum();
    test_colSum();
    test_maxElement();
    test_convergeMatrix();
    test_removeLastRow();
    test_removeLastColumn();
    test_createDiagonalMatrix();
    test_copyMatrix();
    test_getRow();
    test_getColumn();
    test_addRowToMatrix();
    test_addRowOfZeros();
    test_removeRow();
    test_addColumnOfZeros();
    test_removeColumn();
    printf("All matrixUtils tests passed.\n");
    return 0;
}
