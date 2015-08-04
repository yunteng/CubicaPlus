/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <linearalgebra/BLOCK_SPARSE_MATRIX.h>
#if USING_OPENMP
#include <omp.h>
#endif
//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX::BLOCK_SPARSE_MATRIX(int blockRows, int blockCols) :
  _blockRows(blockRows), _blockCols(blockCols)
{
  _blocks = new SpMat[blockRows * blockCols];

  _subRows.resize(_blockRows, 0);
  _subCols.resize(_blockCols, 0);
}

BLOCK_SPARSE_MATRIX::BLOCK_SPARSE_MATRIX() :
  _blocks(NULL), _blockRows(0), _blockCols(0)
{
}

BLOCK_SPARSE_MATRIX::~BLOCK_SPARSE_MATRIX()
{
  if(_blockCols)
    delete[] _blocks;
}

void BLOCK_SPARSE_MATRIX::setBlockDimensions(const vector<int>& rowDims, const vector<int>& colDims)
{
  assert(_blockRows == rowDims.size() && _blockCols == colDims.size());

  _subRows = rowDims;
  _subCols = colDims;
  for(int x = 0; x < _blockRows; x++)
    for(int y = 0; y < _blockCols; y++)
      _blocks[x * _blockCols + y].resize(_subRows[x], _subCols[y]);
}

//////////////////////////////////////////////////////////////////////
// add this matrix to the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::add(const SpMat& matrix, int row, int col)
{
  assert(_subRows[row] == matrix.rows() && _subCols[col] == matrix.cols());

  _blocks[row * _blockCols + col] += matrix;
}

//////////////////////////////////////////////////////////////////////
// subtract this matrix from the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::subtract(const SpMat& matrix, int row, int col)
{
  assert(_subRows[row] == matrix.rows() && _subCols[col] == matrix.cols());

  _blocks[row * _blockCols + col] -= matrix;
}

//////////////////////////////////////////////////////////////////////
// set the block at (row,col) to this matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::equals(const COO_MATRIX& matrix, int row, int col)
{
  assert(_subRows[row] == matrix.rows() && _subCols[col] == matrix.cols());

  matrix.toSpMat(_blocks[row * _blockCols + col]);
  _blocks[row * _blockCols + col].makeCompressed();
}

//////////////////////////////////////////////////////////////////////
// Wipe all entries to zero
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::clear()
{
  int index = 0;
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++, index++)
      _blocks[index].setZero();
}

//////////////////////////////////////////////////////////////////////
// Summed rows across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_SPARSE_MATRIX::totalRows() const
{
  int total = 0;
  for (int x = 0; x < _blockRows; x++)
    total += _subRows[x];
  return total;
}

//////////////////////////////////////////////////////////////////////
// Summed cols across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_SPARSE_MATRIX::totalCols() const
{
  int total = 0;
  for (int x = 0; x < _blockCols; x++)
    total += _subCols[x];
  return total;
}

//////////////////////////////////////////////////////////////////////
// block matrix scalar multiply
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator*=(const Real alpha)
{
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] *= alpha;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::resizeAndWipe(int blockRows, int blockCols)
{
  // clean up old versions
  if(_blockCols)
    delete[] _blocks;
  
  // resize and init
  _blockRows = blockRows;
  _blockCols = blockCols;

  _blocks = new SpMat[blockRows * blockCols];

  for(int x = 0; x < blockRows * blockCols; x++)
    _blocks[x].resize(0, 0);

  _subRows.resize(_blockRows, 0);
  _subCols.resize(_blockCols, 0);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator+=(const BLOCK_SPARSE_MATRIX& A)
{
  assert(_blockRows == A.blockRows());
  assert(_blockCols == A.blockCols());

  for (int x = 0; x < _blockRows; x++){
    for (int y = 0; y < _blockCols; y++){
      add(A(x, y), x, y);
    }
  }

  return *this;
}

BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator-=(const BLOCK_SPARSE_MATRIX& A)
{
  assert(_blockRows == A.blockRows());
  assert(_blockCols == A.blockCols());

  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      subtract(A(x,y), x,y);

  return *this;
}
void BLOCK_SPARSE_MATRIX::mulVec(const VECTOR& v, VECTOR& res) const
{
  assert(totalCols() == v.size());

  res.conservativeResize(totalRows());
  res.setZero();

  vector<int> rowOffsets(_blockRows, 0);
  vector<int> colOffsets(_blockCols, 0);

  for(int x = 1; x < _blockRows; x++)
    rowOffsets[x] = rowOffsets[x - 1] + _subRows[x - 1];
  for(int x = 1; x < _blockCols; x++)
    colOffsets[x] = colOffsets[x - 1] + _subCols[x - 1];

  // int rowOffset = 0;
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int x = 0; x < _blockRows; x++){
    // int colOffset = 0;
    for(int y = 0; y < _blockCols; y++){
      if(_blocks[x * _blockCols + y].nonZeros() == 0)
        continue;
      res.segment(rowOffsets[x], _subRows[x]) += _blocks[x * _blockCols + y] * v.segment(colOffsets[y], _subCols[y]);
      // colOffset += _subCols[y];
    }
    // rowOffset += _subRows[x];
  }
}

void BLOCK_SPARSE_MATRIX::rowMulVec(int row, const VECTOR& v, VECTOR& res) const
{
  assert(totalCols() == v.size());

  res.conservativeResize(_subRows[row]);
  res.setZero();

  vector<int> colOffsets(_blockCols, 0);

  for(int x = 1; x < _blockCols; x++)
    colOffsets[x] = colOffsets[x - 1] + _subCols[x - 1];

  for(int y = 0; y < _blockCols; y++){
    if(_blocks[row * _blockCols + y].nonZeros() == 0)
      continue;
    res += _blocks[row * _blockCols + y] * v.segment(colOffsets[y], _subCols[y]);
  }
}

VECTOR BLOCK_SPARSE_MATRIX::operator*(const VECTOR& v) const
{
  assert(totalCols() == v.size());

  VECTOR res(totalRows());
  res.setZero();

  vector<int> rowOffsets(_blockRows, 0);
  vector<int> colOffsets(_blockCols, 0);

  for(int x = 1; x < _blockRows; x++)
    rowOffsets[x] = rowOffsets[x - 1] + _subRows[x - 1];
  for(int x = 1; x < _blockCols; x++)
    colOffsets[x] = colOffsets[x - 1] + _subCols[x - 1];

#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int x = 0; x < _blockRows; x++){
    for(int y = 0; y < _blockCols; y++){
      if(_blocks[x * _blockCols + y].nonZeros() == 0)
        continue;
      res.segment(rowOffsets[x], _subRows[x]) += _blocks[x * _blockCols + y] * v.segment(colOffsets[y], _subCols[y]);
    }
  }
  return res;
}
void BLOCK_SPARSE_MATRIX::toSpMat(SpMat& spMat)
{
  spMat.resize(totalRows(), totalCols());

  vector<TRIPLET> data;
  int rowOffset = 0;
  for(int x = 0; x < _blockRows; x++){
    int colOffset = 0;
    for(int y = 0; y < _blockCols; y++){
      SpMat& mat = _blocks[x * _blockCols + y];
      for(int k = 0; k < mat.outerSize(); k++){
        for(SpMat::InnerIterator it(mat, k); it; ++it){
          data.push_back(TRIPLET(it.row() + rowOffset, it.col() + colOffset, it.value()));
        }
      }
      colOffset += _subCols[y];
    }
    rowOffset += _subRows[x];
  }
  spMat.setFromTriplets(data.begin(), data.end());
}
