#ifndef BLOCK_COO_MATRIX_H
#define BLOCK_COO_MATRIX_H

#include <SETTINGS.h>
#include <linearalgebra/COO_MATRIX.h>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Block matrix -- essentially a 2D storage class for a bunch of 
// matrices of the same size.
//
// The purpose is to allow matrices to be added and subtracted to
// blocks without having to explicitly check if the blocks exist
// already.
//
// There is no dimension checking, so you better know what you're
// doing.
//////////////////////////////////////////////////////////////////////
class BLOCK_COO_MATRIX {

public:
  BLOCK_COO_MATRIX();

  // create a [blockRows x blockCols] block matrix
  BLOCK_COO_MATRIX(int blockRows, int blockCols);
  
  virtual ~BLOCK_COO_MATRIX();

  void setBlockDimensions(const vector<int>& rowDims, const vector<int>& colDims);

  // set dimensions
  void resizeAndWipe(int blockRows, int blockCols);

  inline COO_MATRIX& operator()(int row, int col) {
    return _blocks[row * _blockCols + col];
  };
  inline const COO_MATRIX& operator()(int row, int col) const { return _blocks[row * _blockCols + col]; };

  BLOCK_COO_MATRIX& operator+=(BLOCK_COO_MATRIX& A);
  BLOCK_COO_MATRIX& operator-=(BLOCK_COO_MATRIX& A);
  BLOCK_COO_MATRIX& operator*=(const Real alpha);
  
  VECTOR operator*(const VECTOR& v) const;

  // add this matrix to the block at (row,col)
  void add(const COO_MATRIX& matrix, int row, int col);

  // subtract this matrix from the block at (row,col)
  void subtract(const COO_MATRIX& matrix, int row, int col);

  // set the block at (row,col) to this matrix
  void equals(const COO_MATRIX& matrix, int row, int col);

  // accessors
  int blockRows() const { return _blockRows; };
  int blockCols() const { return _blockCols; };
  int subRows(int blockRow) const { return _subRows[blockRow]; };
  int subCols(int blockCol) const { return _subCols[blockCol]; };

  // wipe all entries to zero
  void clear();
  
  // summed dims across all blocks
  int totalRows() const;
  int totalCols();

  void toCOOMatrix(COO_MATRIX& cooMat, bool lowerTriangleOnly = false);


protected:
  COO_MATRIX* _blocks;

  // dimensions of the high level 2D block array
  int _blockRows;
  int _blockCols;

  // dimensions of each of the subblocks
  vector<int> _subRows;
  vector<int> _subCols;
};

// BLOCK_VECTOR operator*(BLOCK_COO_MATRIX& A, BLOCK_VECTOR& v);
// VECTOR operator*(const BLOCK_COO_MATRIX& A, const VECTOR& v);
// BLOCK_COO_MATRIX operator+(const BLOCK_COO_MATRIX& A, const BLOCK_COO_MATRIX& B);

#endif
