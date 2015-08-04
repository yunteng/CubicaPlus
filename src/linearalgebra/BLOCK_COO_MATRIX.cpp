/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <linearalgebra/BLOCK_COO_MATRIX.h>
#if USING_OPENMP
#include <omp.h>
#endif
//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_COO_MATRIX::BLOCK_COO_MATRIX(int blockRows, int blockCols) :
  _blockRows(blockRows), _blockCols(blockCols)
{
  _blocks = new COO_MATRIX[blockRows * blockCols];

  _subRows.resize(_blockRows, 0);
  _subCols.resize(_blockCols, 0);
}

BLOCK_COO_MATRIX::BLOCK_COO_MATRIX() :
  _blocks(NULL), _blockRows(0), _blockCols(0)
{
}

BLOCK_COO_MATRIX::~BLOCK_COO_MATRIX()
{
  if(_blockCols)
    delete[] _blocks;
}

void BLOCK_COO_MATRIX::setBlockDimensions(const vector<int>& rowDims, const vector<int>& colDims)
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
void BLOCK_COO_MATRIX::add(const COO_MATRIX& matrix, int row, int col)
{
   assert(_subRows[row] == matrix.rows() && _subCols[col] == matrix.cols());

  _blocks[row * _blockCols + col] += matrix;
}

//////////////////////////////////////////////////////////////////////
// subtract this matrix from the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_COO_MATRIX::subtract(const COO_MATRIX& matrix, int row, int col)
{
  assert(_subRows[row] == matrix.rows() && _subCols[col] == matrix.cols());

  _blocks[row * _blockCols + col] -= matrix;
}

//////////////////////////////////////////////////////////////////////
// set the block at (row,col) to this matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_COO_MATRIX::equals(const COO_MATRIX& matrix, int row, int col)
{
  assert(_subRows[row] == matrix.rows() && _subCols[col] == matrix.cols());

  _blocks[row * _blockCols + col] = matrix;
}

//////////////////////////////////////////////////////////////////////
// Wipe all entries to zero
//////////////////////////////////////////////////////////////////////
void BLOCK_COO_MATRIX::clear()
{
  int index = 0;
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++, index++)
      _blocks[index].clear();
}

//////////////////////////////////////////////////////////////////////
// Summed rows across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_COO_MATRIX::totalRows() const
{
  int total = 0;
  for (int x = 0; x < _blockRows; x++)
    total += _subRows[x];
  return total;
}

//////////////////////////////////////////////////////////////////////
// Summed cols across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_COO_MATRIX::totalCols()
{
  int total = 0;
  for (int x = 0; x < _blockCols; x++)
    total += _subCols[x];
  return total;
}

//////////////////////////////////////////////////////////////////////
// block matrix scalar multiply
//////////////////////////////////////////////////////////////////////
BLOCK_COO_MATRIX& BLOCK_COO_MATRIX::operator*=(const Real alpha)
{
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] *= alpha;

  return *this;
}

void BLOCK_COO_MATRIX::resizeAndWipe(int blockRows, int blockCols)
{
  // clean up old versions
  if(_blockCols)
    delete[] _blocks;
  
  // resize and init
  _blockRows = blockRows;
  _blockCols = blockCols;

  _blocks = new COO_MATRIX[blockRows * blockCols];

  for(int x = 0; x < blockRows * blockCols; x++)
    _blocks[x].resize(0, 0);

  _subRows.resize(_blockRows, 0);
  _subCols.resize(_blockCols, 0);
}

void BLOCK_COO_MATRIX::toCOOMatrix(COO_MATRIX& cooMat, bool lowerTriangleOnly)
{
  cooMat.resize(totalRows(), totalCols());

  vector<TRIPLET>& tuples = cooMat.matrix();

  tuples.clear();

  int rowOffset = 0;

  if(!lowerTriangleOnly){

#if USING_OMP
#pragma omp parallel for schedule(static)
#endif

    for(int x = 0; x < _blockRows; x++)
      for(int y = 0; y < _blockCols; y++){
        _blocks[x * _blockCols + y].order();
        _blocks[x * _blockCols + y].aggregate();
      }

    int index = 0;
    for(int x = 0; x < _blockRows; x++)
    {
      int colOffset = 0;
      for(int y = 0; y < _blockCols; y++, index++)
      {
        vector<TRIPLET>& data = _blocks[index].matrix();
        for(unsigned int z = 0; z < _blocks[index].nnZ(); z++)
          tuples.push_back(TRIPLET(data[z].row() + rowOffset, data[z].col() + colOffset, data[z].value()));
        colOffset += _subCols[y];
      }
      rowOffset += _subRows[x];
    }
  }else{

#if USING_OMP
#pragma omp parallel for schedule(static)
#endif
    for(int x = 0; x < _blockRows; x++)
      for(int y = 0; y <= x; y++){
        _blocks[x * _blockCols + y].order();
        _blocks[x * _blockCols + y].aggregate(x == y);
      }

    for(int x = 0; x < _blockRows; x++)
    {
      int colOffset = 0;
      for(int y = 0; y <= x; y++)
      {
        vector<TRIPLET>& data = _blocks[x * _blockCols + y].matrix();

        for(unsigned int z = 0; z < _blocks[x * _blockCols + y].nnZ(); z++)
          tuples.push_back(TRIPLET(data[z].row() + rowOffset, data[z].col() + colOffset, data[z].value()));

        colOffset += _subCols[y];
      }
      rowOffset += _subRows[x];
    }
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BLOCK_COO_MATRIX& BLOCK_COO_MATRIX::operator+=(BLOCK_COO_MATRIX& A)
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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BLOCK_COO_MATRIX& BLOCK_COO_MATRIX::operator-=(BLOCK_COO_MATRIX& A)
{
  assert(_blockRows == A.blockRows());
  assert(_blockCols == A.blockCols());

  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      subtract(A(x,y), x,y);

  return *this;
}


VECTOR BLOCK_COO_MATRIX::operator*(const VECTOR& v) const
{
  assert(totalRows() == v.size());

  VECTOR res(v.size());
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
      if(_blocks[x * _blockCols + y].nnZ() == 0)
        continue;
      res.segment(rowOffsets[x], _subRows[x]) += _blocks[x * _blockCols + y] * v.segment(colOffsets[y], _subCols[y]);
    }
  }
  return res;
}
