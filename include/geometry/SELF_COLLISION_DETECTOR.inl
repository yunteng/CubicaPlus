#include <geometry/SELF_COLLISION_DETECTOR.h>
#include <util/SIMPLE_PARSER.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#if USING_OPENMP
#include <omp.h>
#endif

template<class BONE>
SELF_COLLISION_DETECTOR<BONE>::SELF_COLLISION_DETECTOR(TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _restSDF(tetMesh->restSDF()),
  _rigger(rigger),
  _skeleton(NULL),
  _surfaceVertexBVH(NULL),
  _tetBVH(NULL),
  _selfCollisionCubatureLoader(NULL)
{
  init(); 
}

template<class BONE>
SELF_COLLISION_DETECTOR<BONE>::~SELF_COLLISION_DETECTOR()
{
  if(_surfaceVertexBVH)
    delete _surfaceVertexBVH;
  if(_tetBVH)
    delete _tetBVH;

  if(_selfCollisionCubatureLoader)
    delete _selfCollisionCubatureLoader;

  for(unsigned int x = 0; x < _partitionedSurfaceVertexBVHs.size(); x++)
    if(_partitionedSurfaceVertexBVHs[x])
      delete _partitionedSurfaceVertexBVHs[x];
  for(unsigned int x = 0; x < _partitionedTetBVHs.size(); x++)
    if(_partitionedTetBVHs[x])
      delete _partitionedTetBVHs[x];
}

template<class BONE>
void SELF_COLLISION_DETECTOR<BONE>::init()
{
  cout << " Initialize Self Collision Detector ..."; flush(cout);
  string sdfFilename = SIMPLE_PARSER::getString("output path", "") + SIMPLE_PARSER::getString("sdf name", "");
  _restSDF.load(sdfFilename);
  _restSDF.initHashSurface(_tetMesh->surfaceFaces());

  _pseudoCCD = SIMPLE_PARSER::getBool("pseudo ccd", false);

  _useLowresTets = SIMPLE_PARSER::getBool("lowres scd", false);
  if(_useLowresTets){
    string lowresMeshName = SIMPLE_PARSER::getString("lowres embedding", "blabla");
    if(!_tetMesh->readLowresEmbeddedMesh(lowresMeshName)){
      cout << " Cannot load lowres embedded mesh at " << lowresMeshName << endl
           << " use high res scd" << endl;
      _useLowresTets = false;
      _useScfCubature = false;
    }else{
      if(!_tetMesh->readLowresEmbeddingMap(lowresMeshName + ".map"))
      {
        _tetMesh->mapToLowresEmbedding();
        _tetMesh->writeLowresEmbeddingMap(lowresMeshName + ".map");
      }
      _useScfCubature = SIMPLE_PARSER::getBool("scf cubature", false);
      if(_useScfCubature){
        _selfCollisionCubatureLoader = new SCF_CUBATURE_LOADER(_tetMesh);
        _useScfCubature = _selfCollisionCubatureLoader->loadAllCubatures(SIMPLE_PARSER::getString("output path", ""));
      }
    }
  }else{
    _useScfCubature = false;
  }

  readSCDList();

  // No Skeleton Information provided
  // do singledomain scd
  if(_rigger == NULL){
    _surfaceVertexBVH = new VERTEX_BVH(_tetMesh);
    _tetBVH = new TET_BVH(_tetMesh, &_restSDF, _useLowresTets);
  }
  else{
    _skeleton = _rigger->skeleton();
    vector<vector<int> > partitionSurfaceVertices;
    vector<vector<int> > partitionTets;

    _rigger->buildSkinningPartition(partitionSurfaceVertices, partitionTets, _useLowresTets);

    _partitionedSurfaceVertexBVHs.resize(_skeleton->totalBones());
    _partitionedTetBVHs.resize(_skeleton->totalBones());

    for(int x = 0; x < _skeleton->totalBones(); x++){

      if(partitionSurfaceVertices[x].size() != 0)
        _partitionedSurfaceVertexBVHs[x] = new VERTEX_BVH(_tetMesh, partitionSurfaceVertices[x]);
      else
        _partitionedSurfaceVertexBVHs[x] = NULL;

      _partitionedTetBVHs[x] = new TET_BVH(_tetMesh, &_restSDF, partitionTets[x], _useLowresTets);
    }
  }

  #if USING_OPENMP
  _selfCollisionPointsCopies.resize(omp_get_max_threads());
  #else
  _selfCollisionPointsCopies.resize(1);
  #endif

  cout << " done. " << endl;
  
}

template<class BONE>
void SELF_COLLISION_DETECTOR<BONE>::readSCDList()
{
  string listFilename = _tetMesh->filename() + ".scd.list";
  FILE* file = fopen(listFilename.c_str(), "r");
  if(file == NULL){
    _selfcheck.resize(_rigger->skeleton()->totalBones(), true);
    _ignoreDomain.resize(_rigger->skeleton()->totalBones(), false);
    _ignoreDomainPair.clear();
    return;
  }
  _selfcheck.resize(_rigger->skeleton()->totalBones(), false);
  _ignoreDomain.resize(_rigger->skeleton()->totalBones(), false);
  _ignoreDomainPair.clear();
  int selfCheckSize = 0;
  fscanf(file, "self_check_domains=%d\n", &selfCheckSize);
  for(int x = 0; x < selfCheckSize; x++){
    int id = 0;
    fscanf(file, "%d\n", &id);
    _selfcheck[id] = true;
    cout << "check self " << id << endl;
  }

  int ignoreDomainSize = 0;
  fscanf(file, "ignore_domains=%d\n", &ignoreDomainSize);
  for(int x = 0; x < ignoreDomainSize; x++){
    int id = 0;
    fscanf(file, "%d\n", &id);
    _ignoreDomain[id] = true;
    // cout << "ignore domain " << id << endl;
  }
  
  int ignoreDomainPairSize = 0;
  fscanf(file, "ignore_pairs=%d\n", &ignoreDomainSize);
  for(int x = 0; x < ignoreDomainSize; x++){
    int id1 = 0;
    int id2 = 0;
    fscanf(file, "%d,%d\n", &id1, &id2);
    _ignoreDomainPair.insert(make_pair(id1, id2));
    _ignoreDomainPair.insert(make_pair(id2, id1));
    cout << "ignore domain pair " << id1 << " " << id2 << endl;
  }
  fclose(file);
}

template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::findNearestSurfacePoint(int vertexPID, VEC3F penetratingRestPosition, SELF_COLLISION_INFO& info, bool checkProjection)
{
  VEC3F nearestRestSurfacePosition;
  if(!_restSDF.nearestSurfacePoint(penetratingRestPosition, nearestRestSurfacePosition)){
    return -1;
  }

  VEC3F* inVertex = _tetMesh->vertex(vertexPID);

  TRIANGLE** surfaceTriangles = NULL;

  int numTris = _restSDF.trianglesAt(nearestRestSurfacePosition, surfaceTriangles);
  // nearest surface point located but failed to find triangle faces
  if(surfaceTriangles == NULL)
    return -1;

  for(int tri = 0; tri < numTris; tri++){
    TRIANGLE*& deformedTriangle = surfaceTriangles[tri];
    int faceID = _tetMesh->surfaceFaceID(deformedTriangle);

    if(checkProjection){
      VEC3F projection = deformedTriangle->projection(*inVertex);
      VEC3F projectionLambda;
      if(!deformedTriangle->baryCenter(projection, projectionLambda))
        continue;
    }

    Real triArea = _tetMesh->restSurfaceFaceArea(surfaceTriangles[tri]);

    VEC3I triPID;
    VEC3F *restVertex[3];
    bool inVertexOnTriangle = false;
    for(int i = 0; i < 3; i++){
      triPID[i] = _tetMesh->vertexID(deformedTriangle->vertices[i]);
      if(triPID[i] == vertexPID){
        inVertexOnTriangle = true;
        break;
      }
      restVertex[i] = _tetMesh->restVertex(triPID[i]);
    }
    if(inVertexOnTriangle)
      continue;

    TRIANGLE restTriangle(*restVertex[0], *restVertex[1], *restVertex[2]);

    VEC3F lambda;
    if(restTriangle.baryCenter(nearestRestSurfacePosition, lambda)){
      Real dist = (penetratingRestPosition - nearestRestSurfacePosition).dot(restTriangle.normal());

      if(dist >= 0)
        continue;

      VEC3F* restPeneratingVertex = _tetMesh->restVertex(vertexPID);
      VEC3F triangleCenter = restTriangle.center();
      Real edge = (*restVertex[0] - *restVertex[1]).norm();
      edge = max(edge, (*restVertex[0] - *restVertex[2]).norm());
      edge = max(edge, (*restVertex[1] - *restVertex[2]).norm());

      Real restDist = (*restPeneratingVertex - nearestRestSurfacePosition).norm();
      if(restDist <= edge * 1.8)
        continue;

      info.vertexID = vertexPID;
      info.triangleVertexIDs = triPID;
      info.faceID = faceID;
      info.baryCenter = lambda;
     
      Real area = (_tetMesh->restSurfaceVertexArea(inVertex) + triArea) / 4.0;

      info.avgArea = area;
      
      return (penetratingRestPosition - nearestRestSurfacePosition).norm();
    }
  }
  return -1;
}
template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::singleDomainVertexVsTetSCD()
{
  Real penetrationDepth = 0;
  int cnt = 0;

  _selfCollisionPoints.clear();

  if(_useLowresTets)
    _tetMesh->updateLowresEmbedding();

  _surfaceVertexBVH->refit();
  _tetBVH->refit();

  vector<pair<int, VEC3F> > collidingNodes;
  _tetBVH->collide(_surfaceVertexBVH, collidingNodes);

  set<int> inCollisionVertices;
  
  for(unsigned int i = 0; i < collidingNodes.size(); i++){
    int vertexPID = collidingNodes[i].first;
    if(inCollisionVertices.find(vertexPID) != inCollisionVertices.end())
        continue;

    SELF_COLLISION_INFO info;
    Real dist = findNearestSurfacePoint(vertexPID, collidingNodes[i].second, info);
    if(dist > 0){
      inCollisionVertices.insert(vertexPID);
      _selfCollisionPoints.push_back(info);
      cnt++;
      penetrationDepth += dist;
    }
  }

  if(cnt == 0)
    cout << "vertexVsTetSCD: self collision free" << endl;
  else{
    cout << "vertexVsTetSCD: " << cnt << " collision pairs " << " avg penetration depth " << penetrationDepth / cnt << endl;
  }

  if(cnt > 0)
    return penetrationDepth;
  return -1;
}

template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::multiDomainVertexVsTetSCD()
{
  int cnt = 0;
  Real penetrationDepth = 0;

  _selfCollisionPoints.clear();

  if(_useLowresTets)
    _tetMesh->updateLowresEmbedding();

  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for(int x = 0; x < _skeleton->totalBones(); x++){
    if(_partitionedSurfaceVertexBVHs[x])
      _partitionedSurfaceVertexBVHs[x]->refit();
  }
  TIMING_BREAKDOWN::toc("Refit Surface BVH");

  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for(int x = 0; x < _skeleton->totalBones(); x++){
    _partitionedTetBVHs[x]->refit();
  }
  TIMING_BREAKDOWN::toc("Refit VOLUME BVH");

#if USING_OPENMP
  TIMING_BREAKDOWN::tic();
  for(unsigned int x = 0; x < _selfCollisionPointsCopies.size(); x++)
    _selfCollisionPointsCopies[x].clear();
  
  int shift = SIMPLE_PARSER::getBool("self scd", true) ? 0 : 1;

  for(; shift < _skeleton->totalBones(); shift++){
    Real localPenetrationDepth = 0;
    int localCnt = 0;
    #pragma omp parallel
    {
      const int id = omp_get_thread_num();

      #pragma omp for schedule(dynamic) reduction(+ : localPenetrationDepth, localCnt)
      for(int x = 0; x < _skeleton->totalBones(); x++){
        int y = (x + shift) % _skeleton->totalBones();

        if(_partitionedSurfaceVertexBVHs[y] == NULL)
          continue;

        if(x == y && !_selfcheck[x])
          continue;

        if(_ignoreDomain[x] || _ignoreDomain[y])
          continue;

        if(_ignoreDomainPair.find(make_pair(x, y)) != _ignoreDomainPair.end())
          continue;

        vector<pair<int, VEC3F> > collidingNodes;
        _partitionedTetBVHs[x]->collide(_partitionedSurfaceVertexBVHs[y], collidingNodes);

        for(unsigned int i = 0; i < collidingNodes.size(); i++){
          int vertexPID = collidingNodes[i].first;
          SELF_COLLISION_INFO info;
          Real dist = findNearestSurfacePoint(vertexPID, collidingNodes[i].second, info, x == y);
          if(dist > 0){
            info.vertexPartition = y;
            info.trianglePartition = x;
            info.cubatureWeight = 1.0;
            info.isCubaturePair = false;

            _selfCollisionPointsCopies[id].push_back(info);
            localCnt++;
            localPenetrationDepth += dist;
          }
        }
      }
    }
    cnt += localCnt;
    penetrationDepth += localPenetrationDepth;
  }
  set<int> inCollisionVertices;
  for(unsigned int x = 0; x < _selfCollisionPointsCopies.size(); x++){
    for(unsigned int y = 0; y < _selfCollisionPointsCopies[x].size(); y++){
      if(inCollisionVertices.find(_selfCollisionPointsCopies[x][y].vertexID) == inCollisionVertices.end()){
        inCollisionVertices.insert(_selfCollisionPointsCopies[x][y].vertexID);
        SELF_COLLISION_INFO& info = _selfCollisionPointsCopies[x][y];

        _selfCollisionPoints.push_back(info);
      }
    }
  }

  TIMING_BREAKDOWN::toc("Check Collisions");
#else
  TIMING_BREAKDOWN::tic();
  int shift = SIMPLE_PARSER::getBool("self scd", true) ? 0 : 1;
  set<int> inCollisionVertices;
  for(; shift < _skeleton->totalBones(); shift++){
    for(int x = 0; x < _skeleton->totalBones(); x++){
      int y = (x + shift) % _skeleton->totalBones();

      if(x == y && !_selfcheck[x])
        continue;

      if(_ignoreDomain[x] || _ignoreDomain[y])
        continue;

      if(_ignoreDomainPair.find(make_pair(x, y)) != _ignoreDomainPair.end())
        continue;

      vector<pair<int, VEC3F> > collidingNodes;
      _partitionedTetBVHs[x]->collide(_partitionedSurfaceVertexBVHs[y], collidingNodes);

      for(unsigned int i = 0; i < collidingNodes.size(); i++){
        int vertexPID = collidingNodes[i].first;
        if(inCollisionVertices.find(vertexPID) != inCollisionVertices.end())
            continue;

        SELF_COLLISION_INFO info;
        Real dist = findNearestSurfacePoint(vertexPID, collidingNodes[i].second, info);
        if(dist > 0){
            inCollisionVertices.insert(vertexPID);
            info.isCubaturePair = false;
            _selfCollisionPoints.push_back(info);  
            cnt++;
            penetrationDepth += dist;
        }
      }
    }
  }
  TIMING_BREAKDOWN::toc("Check Collisions");
#endif

  if(cnt == 0)
    cout << "vertexVsTetSCD: self collision free" << endl;
  else{
    cout << "vertexVsTetSCD: " << cnt << " collision pairs " << " avg penetration depth " << penetrationDepth / cnt << endl;
  }

  if(cnt > 0)
    return penetrationDepth;
  return -1;
}

template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::lowresCubatureSelfCollisionTest()
{
  cout << "============Low res cubature self collision test============" << endl;

  int cnt = 0;
  Real penetrationDepth = 0;

  _selfCollisionPoints.clear();
  
  _tetMesh->updateLowresEmbedding();

  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for(int x = 0; x < _skeleton->totalBones(); x++){
    if(_partitionedSurfaceVertexBVHs[x])
      _partitionedSurfaceVertexBVHs[x]->refit();
  }
  TIMING_BREAKDOWN::toc("Refit Surface BVH");

  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for(int x = 0; x < _skeleton->totalBones(); x++){
    _partitionedTetBVHs[x]->refit();
  }
  TIMING_BREAKDOWN::toc("Refit VOLUME BVH");

  for(int x = 0; x < _skeleton->totalBones(); x++)
    for(int y = x + 1; y < _skeleton->totalBones(); y++){
      VEC3F t;
      QUATERNION r;
      _skeleton->computeRelativeTransform(x, y, t, r);
      _relativeTransforms[make_pair(x, y)] = make_pair(t, r);
      _relativeTransforms[make_pair(y, x)] = make_pair(t, r);
    }

  TIMING_BREAKDOWN::tic();
  for(unsigned int x = 0; x < _selfCollisionPointsCopies.size(); x++)
    _selfCollisionPointsCopies[x].clear();

  int shift = 1;
  for(; shift < _skeleton->totalBones(); shift++){
    Real localPenetrationDepth = 0;
    int localCnt = 0;

    #if USING_OPENMP
    #pragma omp parallel
    #endif
    {
      #if USING_OPENMP
      const int id = omp_get_thread_num();
      #else
      const int id = 0;
      #endif

      #if USING_OPENMP
      #pragma omp for schedule(static) reduction(+ : localPenetrationDepth, localCnt)
      #endif
      for(int x = 0; x < _skeleton->totalBones(); x++){
        int y = (x + shift) % _skeleton->totalBones();

        if(_partitionedSurfaceVertexBVHs[y] == NULL)
          continue;

        if(_ignoreDomain[x] || _ignoreDomain[y])
          continue;

        if(_ignoreDomainPair.find(make_pair(x, y)) != _ignoreDomainPair.end())
          continue;

        pair<int, int> partitionPair(x, y);

        vector<pair<int, Real> > leftVertexIDs, rightVertexIDs;

        vector<VEC3F> rightNodes;
        vector<VEC3F> leftNodes;

        // if the training pose space has captured this potential collision pair, just do fast patch-volume intersection test
        pair<VEC3F, QUATERNION>& relativeTransform = _relativeTransforms[partitionPair];
        const VEC3F& relativeTranslation = relativeTransform.first;
        const QUATERNION& relativeRotation = relativeTransform.second;

        bool cacheFound = false;
        #if USING_OPENMP
        #pragma omp critical
        #endif
        {
          cacheFound = _selfCollisionCubatureLoader->getCubatureSurfaceVertices(partitionPair, relativeTranslation, relativeRotation, leftVertexIDs, rightVertexIDs);
        }

        if(cacheFound){
          if(_tetMesh->isNeighbors(x, y))
            for(unsigned int z = 0; z < leftVertexIDs.size(); z++)
              leftNodes.push_back(*(_tetMesh->vertex(leftVertexIDs[z].first)));

          for(unsigned int z = 0; z < rightVertexIDs.size(); z++)
            rightNodes.push_back(*(_tetMesh->vertex(rightVertexIDs[z].first)));
        }

        if(!leftNodes.empty()){
          vector<pair<int, VEC3F> > collidingNodes;
          _partitionedTetBVHs[x]->collide(leftNodes, collidingNodes);

          for(unsigned int i = 0; i < collidingNodes.size(); i++){
            int vertexPID = leftVertexIDs[collidingNodes[i].first].first;
            Real cubatureWeight = leftVertexIDs[collidingNodes[i].first].second;
            
            SELF_COLLISION_INFO info;
            Real dist = findNearestSurfacePoint(vertexPID, collidingNodes[i].second, info, true);
            if(dist > 0){
              info.cubatureWeight = cubatureWeight;
              info.vertexPartition = y;
              info.trianglePartition = x;
              info.isCubaturePair = true;

              _selfCollisionPointsCopies[id].push_back(info);
              localCnt++;
              localPenetrationDepth += dist;
            }
          }
        }

        if(!rightNodes.empty()){
          vector<pair<int, VEC3F> > collidingNodes;
          _partitionedTetBVHs[x]->collide(rightNodes, collidingNodes);

          for(unsigned int i = 0; i < collidingNodes.size(); i++){
            int vertexPID = rightVertexIDs[collidingNodes[i].first].first;
            Real cubatureWeight = rightVertexIDs[collidingNodes[i].first].second;

            SELF_COLLISION_INFO info;
            info.vertexPartition = y;
            info.trianglePartition = x;

            Real dist = findNearestSurfacePoint(vertexPID, collidingNodes[i].second, info, false);
            if(dist > 0){

              info.cubatureWeight = cubatureWeight;
              info.isCubaturePair = true;

              _selfCollisionPointsCopies[id].push_back(info); 
              localCnt++;
              localPenetrationDepth += dist;
            }
          }
        }
        if(!cacheFound){
          vector<pair<int, VEC3F> > collidingNodes;
          _partitionedTetBVHs[x]->collide(_partitionedSurfaceVertexBVHs[y], collidingNodes);
          for(unsigned int i = 0; i < collidingNodes.size(); i++){
            int vertexPID = collidingNodes[i].first;
            SELF_COLLISION_INFO info;
            info.vertexPartition = y;
            info.trianglePartition = x;
            info.cubatureWeight = 1.0;

            Real dist = findNearestSurfacePoint(vertexPID, collidingNodes[i].second, info, x == y);
            if(dist > 0){
              info.isCubaturePair = false;

              _selfCollisionPointsCopies[id].push_back(info);
              localCnt++;
              localPenetrationDepth += dist;
            }
          }
        }       
      }
      cnt += localCnt;
      penetrationDepth += localPenetrationDepth;
    }
  }

  set<int> inCollisionVertices;
  for(unsigned int x = 0; x < _selfCollisionPointsCopies.size(); x++){
    for(unsigned int y = 0; y < _selfCollisionPointsCopies[x].size(); y++){
      if(inCollisionVertices.find(_selfCollisionPointsCopies[x][y].vertexID) == inCollisionVertices.end()){
        inCollisionVertices.insert(_selfCollisionPointsCopies[x][y].vertexID);
        SELF_COLLISION_INFO& info = _selfCollisionPointsCopies[x][y];

        _selfCollisionPoints.push_back(info);
      }
    }
  }

  TIMING_BREAKDOWN::toc("Check Collisions");

  if(cnt == 0)
    cout << "Cubature SCD: self collision free" << endl;
  else{
    cout << "Cubature SCD: " << _selfCollisionPoints.size() << " collision pairs " << " avg penetration depth " << penetrationDepth / cnt << endl;
  }
  if(cnt > 0)
    return penetrationDepth;
  
  return -1;
}

template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::vertexVsTetSCD()
{
  static int frameCnt = 25;

  if(_rigger == NULL)
    return singleDomainVertexVsTetSCD();
  else{
    if(_pseudoCCD){
      VECTOR& originalDisp = _tetMesh->x();
      VECTOR& skinningDisp = _rigger->skinningDisp();

      Real w = SIMPLE_PARSER::getFloat("ccd weight", 0.1);

      VECTOR diff = skinningDisp - originalDisp;
      vector<VEC3F>& vertices = _tetMesh->vertices();

      for(unsigned int x = 0; x < diff.size() / 3; x++){
        vertices[x] += w * diff.segment<3>(x * 3); 
      }

    }
    if(!_useScfCubature)
      return multiDomainVertexVsTetSCD();
    else
      return lowresCubatureSelfCollisionTest();
    if(_pseudoCCD)
      _tetMesh->updateFullMesh();
  }
  frameCnt++;
  if(frameCnt > 105)
    frameCnt = 34;
}

template<class BONE>
void SELF_COLLISION_DETECTOR<BONE>::drawSelfCollisionPoints()
{
  glDisable(GL_DEPTH_TEST);

  for(unsigned int x = 0; x < _selfCollisionPoints.size(); x++){
    SELF_COLLISION_INFO& info = _selfCollisionPoints[x];

    VEC3F vertex = *_tetMesh->vertex(info.vertexID);

    VEC3I& triangleIndex = info.triangleVertexIDs;
    VEC3F& lambda = info.baryCenter;

    VEC3F surfacePosition;
    surfacePosition.setZero();
    for(int x = 0; x < 3; x++){
      surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndex[x]));
    }
    
    glPointSize(5.0f);

    glBegin(GL_POINTS);
    glColor4f(0.0, 0.0, 1.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glColor4f(1.0, 0.0, 0.0, 1.0);
      glVertex3f(surfacePosition[0], surfacePosition[1], surfacePosition[2]);
    glEnd();

    glLineWidth(2.0f);
    glColor4f(10.0, 10.0, 0.0, 1.0);
    glBegin(GL_LINES);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(surfacePosition[0], surfacePosition[1], surfacePosition[2]);
    glEnd();
  }
  glEnable(GL_DEPTH_TEST);
}

template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::collide(vector<VEC3F>& collidingNodes)
{
  Real penetrationDepth = 0;
  int cnt = 0;

  _externalCollisionPoints.clear();
  
  for(unsigned int i = 0; i < collidingNodes.size(); i++){

    EX_COLLISION_INFO info;
    Real dist = findNearestSurfacePoint(collidingNodes[i], info);
    if(dist > 0){
      info.penetratingPosition = collidingNodes[i];

      _externalCollisionPoints.push_back(info);
      cnt++;
      penetrationDepth += dist;
    }
  }

  if(cnt != 0){
    cout << "collide: " << cnt << " collision pairs " << " avg penetration depth " << penetrationDepth / cnt << endl;
  }

  if(cnt > 0)
    return penetrationDepth;
  return -1;
}

template<class BONE>
Real SELF_COLLISION_DETECTOR<BONE>::findNearestSurfacePoint(const VEC3F& penetratingRestPosition, EX_COLLISION_INFO& info)
{
  VEC3F nearestRestSurfacePosition;
  if(!_restSDF.nearestSurfacePoint(penetratingRestPosition, nearestRestSurfacePosition)){
    return -1;
  }

  TRIANGLE** surfaceTriangles = NULL;

  int numTris = _restSDF.trianglesAt(nearestRestSurfacePosition, surfaceTriangles);
  // nearest surface point located but failed to find triangle faces
  if(surfaceTriangles == NULL)
    return -1;

  for(int tri = 0; tri < numTris; tri++){
    TRIANGLE*& deformedTriangle = surfaceTriangles[tri];
    int faceID = _tetMesh->surfaceFaceID(deformedTriangle);

    VEC3I triPID;
    VEC3F *restVertex[3];

    for(int i = 0; i < 3; i++){
      triPID[i] = _tetMesh->vertexID(deformedTriangle->vertices[i]);
      restVertex[i] = _tetMesh->restVertex(triPID[i]);
    }

    TRIANGLE restTriangle(*restVertex[0], *restVertex[1], *restVertex[2]);
   
    VEC3F lambda;
    if(restTriangle.baryCenter(nearestRestSurfacePosition, lambda)){
      Real dist = (penetratingRestPosition - nearestRestSurfacePosition).dot(restTriangle.normal());

      if(dist >= 0)
        continue;

      info.triangleVertexIDs = triPID;
      info.baryCenter = lambda;
      info.avgArea = _tetMesh->restSurfaceFaceArea(surfaceTriangles[tri]) / 3.0;
     
      return (penetratingRestPosition - nearestRestSurfacePosition).norm();
    }
  }
  return -1;
}
template<class BONE>
void SELF_COLLISION_DETECTOR<BONE>::drawExternalCollisionPoints()
{
  glDisable(GL_DEPTH_TEST);

  for(unsigned int x = 0; x < _externalCollisionPoints.size(); x++){
    EX_COLLISION_INFO& info = _externalCollisionPoints[x];

    VEC3F vertex = info.penetratingPosition;

    VEC3I& triangleIndex = info.triangleVertexIDs;
    VEC3F& lambda = info.baryCenter;

    VEC3F surfacePosition;
    surfacePosition.setZero();
    for(int x = 0; x < 3; x++){
      surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndex[x]));
    }
    
    glPointSize(5.0f);

    glBegin(GL_POINTS);
    glColor4f(0.0, 0.0, 1.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glColor4f(1.0, 0.0, 0.0, 1.0);
      glVertex3f(surfacePosition[0], surfacePosition[1], surfacePosition[2]);
    glEnd();

    glLineWidth(2.0f);
    glColor4f(10.0, 10.0, 0.0, 1.0);
    glBegin(GL_LINES);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(surfacePosition[0], surfacePosition[1], surfacePosition[2]);
    glEnd();
  }
  glEnable(GL_DEPTH_TEST);
}
