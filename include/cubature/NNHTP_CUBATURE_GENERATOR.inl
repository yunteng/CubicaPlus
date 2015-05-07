#include <util/NNLS.h>
template<class APP>
NNHTP_CUBATURE_GENERATOR<APP>::NNHTP_CUBATURE_GENERATOR(APP* app):
  _app(app),
  _trainingTwister(654321)
{

}
template<class APP>
NNHTP_CUBATURE_GENERATOR<APP>::~NNHTP_CUBATURE_GENERATOR()
{
  
}
template<class APP>
void NNHTP_CUBATURE_GENERATOR<APP>::generateCubatures()
{
  if(_app->numberOfSamples() == 0){
    cout << " NNHTP cubature generator: No samples to train!!!" << endl;
    return;
  }
  
  // size of a single force vector
  int forceSize = _app->sampleDimension();
  
  // total number of force samples
  int totalSamples = _app->numberOfSamples();

  // number of rows in a column where all the samples are flattened
  // into a big vector
  int totalRows = forceSize * totalSamples;


  // stack force samples into one big 'b' vector
  VECTOR b(totalRows);
  int index = 0;
  for (int x = 0; x < totalSamples; x++)
    for (int y = 0; y < forceSize; y++, index++)
      b(index) = _app->samples()[x](y);

  _app->clearSamples();

  int totalCandidates = _app->totalCandidates();
  int maxKeyPoints = _app->maxKeyPoints();
  Real errorTolerance = _app->errorTolerance();

  cout << " total candidates " << totalCandidates << endl
       << " max key points " << maxKeyPoints << endl
       << " error tolerance " << errorTolerance << endl;

  VECTOR gradient(totalCandidates);
  gradient.setZero();

  VECTOR gradientNoW(totalCandidates);
  gradientNoW.setZero();

  VECTOR gradientS(totalCandidates);
  gradientS.setZero();

  VECTOR w(totalCandidates);
  w.setZero();

  VECTOR oldW(totalCandidates);
  oldW.setZero();

  vector<int> wIndex;

  vector<int>& initialKeyPoints = _app->initialKeyPoints();
  vector<Real>& initialWeights = _app->initialWeights();

  // if no precomputed cubatures are found, just random initialize the cubatures and set all weights to 1
  if(initialKeyPoints.empty()){
    wIndex = NNHTP_randomPickCandidates(wIndex, maxKeyPoints, totalCandidates, false);
    for(unsigned int x = 0; x < wIndex.size(); x++)
      w[wIndex[x]] = 1.0;
    cout << " Random Initial Guess, picked " << wIndex.size() << " points" << endl;
  }else{
    cout << " Initial Guess found" << endl;

    _keyPointIDs = initialKeyPoints;
    _keyWeights = initialWeights;
    verifyCubaturePoints(b);

    wIndex = initialKeyPoints;
    for(unsigned int x = 0; x < wIndex.size(); x++)
      w[wIndex[x]] = initialWeights[x];
  }
  

  int maxIteration = _app->maxIteration();

  NNLS_SOLVER nnls(totalRows, maxKeyPoints);
  nnls.maxIter() = 10;

  double* bNNLS = new double[totalRows];
  double* weightsNNLS = new double[maxKeyPoints];
  memset(weightsNNLS, 0, maxKeyPoints * sizeof(double));

  // "A" matrix to send to the the NNLS solver
  double* A = new double[totalRows * maxKeyPoints];
  double bNorm = b.norm();
  double rNorm = bNorm;
  Real relativeError = 1.0;

  for(int i = 0; i < maxIteration; i++){

    oldW = w;

    // lazy evaluate the gradient
    vector<int> candidateSubset = NNHTP_randomPickCandidates(wIndex, maxKeyPoints * 5, totalCandidates, true);
    map<int, VECTOR> sampleColumns;

    // gradient = 2A^T(Aw-b)
    for(unsigned int x = 0; x < candidateSubset.size(); x++){
      sampleColumns[candidateSubset[x]] = _app->getCandidateQuantitiy(candidateSubset[x]);
    }
    VECTOR Aw(totalRows);
    Aw.setZero();
    for(unsigned int x = 0; x < wIndex.size(); x++)
      Aw += (w[wIndex[x]] * sampleColumns[wIndex[x]]);

    Aw -= b;

    gradient.setZero();
// #pragma omp parallel for
    for(unsigned int x = 0; x < candidateSubset.size(); x++){
      gradient[candidateSubset[x]] = 2 * (sampleColumns[candidateSubset[x]].dot(Aw));
    }

    gradientS = gradient;
    VECTOR Aw_s(totalRows);
    Aw_s.setZero();
    for(unsigned int x = 0; x < candidateSubset.size(); x++){
      if(gradientS[candidateSubset[x]] > 0)
        Aw_s += (gradientS[candidateSubset[x]] * sampleColumns[candidateSubset[x]]);
      else
        gradientS[candidateSubset[x]] = 0;
    }

    Real stepSize = gradientS.squaredNorm() / Aw_s.squaredNorm();
    if(std::isnan(stepSize))
      continue;
    // cout << "step size " << stepSize << endl;
    w -= (stepSize * gradient);

    vector<int> newCandidates = NNHTP_project(w, maxKeyPoints);
    
    for(unsigned int x = 0; x < newCandidates.size(); x++){
      if(sampleColumns.find(newCandidates[x]) == sampleColumns.end()){
        
        sampleColumns[newCandidates[x]] = _app->getCandidateQuantitiy(newCandidates[x]);
      }
      VECTOR& column = sampleColumns[newCandidates[x]];
      for(unsigned int y = 0; y < totalRows; y++){
        A[x * totalRows + y] = column(y);
      }
    }
    for (int x = 0; x < totalRows; x++)
      bNNLS[x] = b(x);

    bool converged = nnls.solve(A, newCandidates.size(), bNNLS, weightsNNLS, rNorm);

    relativeError = rNorm / bNorm;
    
    cout << "iteration " << i << " relative error " << relativeError << endl;

    w.setZero();
    wIndex.clear();
    for(int x = 0; x < newCandidates.size(); x++){
      if(weightsNNLS[x] <= 1e-12)
        continue;

      int idx = newCandidates[x];
      w[idx] = weightsNNLS[x];
      wIndex.push_back(idx);
    }
    if(relativeError <= errorTolerance)
      break;
    
    if((w - oldW).squaredNorm() < 1e-6){
      cout << "picked same samples in two consecutive steps!!! stop" << endl;
      break;
    }
    // save some partial result so we don't have to wait forever
    if(i > 0 && i % 5 == 0){
      _keyPointIDs.clear();
      _keyWeights.clear();
      for(unsigned int x = 0; x < wIndex.size(); x++){
        if(w[wIndex[x]] < 1e-12)
          continue;
        _keyPointIDs.push_back(wIndex[x]);
        _keyWeights.push_back(w[wIndex[x]]);
      }
      cout << " wrtie partial result" << endl;
      writeCubatures(_app->cubatureFilename() + ".iteration." + IO::intToString(i));
    }
  }

  _keyPointIDs.clear();
  _keyWeights.clear();

  for(unsigned int x = 0; x < wIndex.size(); x++){
    if(w[wIndex[x]] < 1e-12)
      continue;
    _keyPointIDs.push_back(wIndex[x]);
    _keyWeights.push_back(w[wIndex[x]]);
  }
  cout << "generated " << _keyPointIDs.size() << " cubatures with relative error " << relativeError << endl;

  verifyCubaturePoints(b);

  delete[] bNNLS;
  delete[] weightsNNLS;
  delete[] A;

  writeCubatures(_app->cubatureFilename());
}

template<class APP>
vector<int> NNHTP_CUBATURE_GENERATOR<APP>::NNHTP_project(VECTOR& w, int toKeep)
{
  priority_queue<pair<int, Real>, vector<pair<int, Real> >, orderByValueSmallerThan> candidateHeap;
  for(int x = 0; x < w.size(); x++){
    if(w[x] > 0)
      candidateHeap.push(make_pair(x, w[x]));
  }
  w.setZero();
  vector<int> indices;
  toKeep = candidateHeap.size() < toKeep ? candidateHeap.size() : toKeep;
  for(int x = 0; x < toKeep; x++){
    pair<int, Real> candidate = candidateHeap.top();

    indices.push_back(candidate.first);
    w[candidate.first] = candidate.second;
    candidateHeap.pop();
  }
  return indices;
}

template<class APP>
vector<int> NNHTP_CUBATURE_GENERATOR<APP>::NNHTP_randomPickCandidates(vector<int>& excludes, int candidatesPerTry, int totalCandidates, bool addExcludes)
{
  vector<bool> alreadyused(totalCandidates, false);
  for(unsigned int x = 0; x < excludes.size(); x++)
    alreadyused[excludes[x]] = true;

  MERSENNETWISTER& twister = _trainingTwister;

  vector<int> candidates;

  if(addExcludes)
    candidates = excludes;

  candidatesPerTry = candidatesPerTry < (totalCandidates - candidates.size()) ? candidatesPerTry : (totalCandidates - candidates.size());

  while(candidates.size() < candidatesPerTry){
    int index = twister.randInt(totalCandidates - 1);
    if(!alreadyused[index]){
      candidates.push_back(index);
      alreadyused[index] = true;
    }
  }
  return candidates;
}

template<class APP>
void NNHTP_CUBATURE_GENERATOR<APP>::verifyCubaturePoints(const VECTOR& b)
{ 
  VECTOR residual = b;

  int cnt = _keyPointIDs.size();

  priority_queue<pair<int, Real>, vector<pair<int, Real> >, orderByValueSmallerThan> candidateHeap;

  for(unsigned int y = 0; y < _keyPointIDs.size(); y++){
    candidateHeap.push(make_pair(_keyPointIDs[y], _keyWeights[y]));
  }
  _keyPointIDs.clear();
  _keyWeights.clear();

  while(!candidateHeap.empty()){
    pair<int, Real> candidate = candidateHeap.top();
    residual -= candidate.second* _app->getCandidateQuantitiy(candidate.first);
    _keyPointIDs.push_back(candidate.first);
    _keyWeights.push_back(candidate.second);

    if(residual.norm() / b.norm() < _app->errorTolerance())
      break;

    candidateHeap.pop();
  }

  // for(unsigned int y = 0; y < _keyPointIDs.size(); y++){
    // residual -= _keyWeights[y] * _app->getCandidateQuantitiy(_keyPointIDs[y]);
    // if(residual.norm() / b.norm() < _app->errorTolerance())
      // break;
  // }

  cout << "verify cubature points: relative error " << residual.norm() / b.norm() << endl;
  cout << "using " << _keyPointIDs.size() << " out of " << cnt << " cubature points" << endl;
}

template<class APP>
void NNHTP_CUBATURE_GENERATOR<APP>::writeCubatures(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "wb");
  if(file == NULL){
    cout << " cannot open " << filename << " to write!!" << endl;
    return;
  }
  int size = _keyPointIDs.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  for(unsigned int x = 0; x < _keyPointIDs.size(); x++){
    int id = _app->getTrueID(_keyPointIDs[x]);
    double weight = _keyWeights[x];
    fwrite((void*)&id, sizeof(int), 1, file);
    fwrite((void*)&weight, sizeof(double), 1, file);
  }

  fclose(file);
}
