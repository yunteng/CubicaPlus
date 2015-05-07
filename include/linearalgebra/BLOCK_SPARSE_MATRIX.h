#ifndef BLOCK_SPARSE_MATRIX_H
#define BLOCK_SPARSE_MATRIX_H

#include <SETTINGS.h>
#include <iostream>
#include <linearalgebra/COO_MATRIX.h>

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
class BLOCK_SPARSE_MATRIX {

public:
  BLOCK_SPARSE_MATRIX();

  // create a [blockRows x blockCols] block matrix
  BLOCK_SPARSE_MATRIX(int blockRows, int blockCols);
  
  virtual ~BLOCK_SPARSE_MATRIX();

  void setBlockDimensions(const vector<int>& rowDims, const vector<int>& colDims);

  // set dimensions
  void resizeAndWipe(int blockRows, int blockCols);
  
  inline SpMat& operator()(int row, int col) {
    return _blocks[row * _blockCols + col];
  };

  inline const SpMat& operator()(int row, int col) const { return _blocks[row * _blockCols + col]; };

  BLOCK_SPARSE_MATRIX& operator+=(const BLOCK_SPARSE_MATRIX& A);
  BLOCK_SPARSE_MATRIX& operator-=(const BLOCK_SPARSE_MATRIX& A);
  BLOCK_SPARSE_MATRIX& operator*=(const Real alpha);
  
  VECTOR operator*(const VECTOR& v) const;
  void mulVec(const VECTOR& v, VECTOR& res) const;
  void rowMulVec(int row, const VECTOR& v, VECTOR& res) const;

  // add this matrix to the block at (row,col)
  void add(const SpMat& matrix, int row, int col);

  // subtract this matrix from the block at (row,col)
  void subtract(const SpMat& matrix, int row, int col);

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
  int totalCols() const;

  void toSpMat(SpMat& cooMat);


protected:
  SpMat* _blocks;

  // dimensions of the high level 2D block array
  int _blockRows;
  int _blockCols;

  // dimensions of each of the subblocks
  vector<int> _subRows;
  vector<int> _subCols;
};

// BLOCK_VECTOR operator*(BLOCK_SPARSE_MATRIX& A, BLOCK_VECTOR& v);
// VECTOR operator*(const BLOCK_SPARSE_MATRIX& A, const VECTOR& v);
// BLOCK_SPARSE_MATRIX operator+(const BLOCK_SPARSE_MATRIX& A, const BLOCK_SPARSE_MATRIX& B);

#endif
