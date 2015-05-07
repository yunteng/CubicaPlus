#include <geometry/TET_MESH.h>
#include <util/RGB_HSV.h>
void TET_MESH::drawSurfaceFaces()
{
  // glColor3f(0.498, 0.784, 0.898);
  for (unsigned int x = 0 ; x < _surfaceFaces.size(); x++){
    _surfaceFaces[x].draw();
  }
}
void TET_MESH::drawPartition(int partition)
{
  for(unsigned int x = 0; x < _surfaceFaces.size(); x++){
    TRIANGLE& face = _surfaceFaces[x];
    int indices[] = {_vertexID[face.vertices[0]],
                     _vertexID[face.vertices[1]],
                     _vertexID[face.vertices[2]]};
    if(_partitionIDs[indices[0]].first == partition ||
       _partitionIDs[indices[1]].first == partition ||
       _partitionIDs[indices[2]].first == partition)
      face.draw();
  }
}
void TET_MESH::drawInterfaceVertices()
{
  glDisable(GL_DEPTH_TEST);
  glPointSize(5.0f);
  glColor3f(1.0, 1.0, 0.0);
  glBegin(GL_POINTS);
  for(unsigned int x = 0; x < _surfaceVertexSize; x++){
    if(_vertexNumberOfCopies[x] > 1){
      glVertex3f(_vertices[x][0], _vertices[x][1], _vertices[x][2]);
    }
  }
  glEnd();
  glEnable(GL_DEPTH_TEST);
}

void TET_MESH::drawFullsimVertices()
{
  glDisable(GL_DEPTH_TEST);
  glPointSize(5.0f);
  glColor3f(1.0, 1.0, 0.0);
  glBegin(GL_POINTS);
  for(int x = 0; x < _totalPartitions; x++){
    for(unsigned int y = 0; y < _partitionedVertices[x].size(); y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];

      if(index * 3 < _partitionedFullsimDofs[x]){
        int vid = _partitionedVertices[x][y];
        glVertex3f(_vertices[vid][0], _vertices[vid][1], _vertices[vid][2]);
      }
    }
  }
  glEnd();
  glEnable(GL_DEPTH_TEST);
}

void TET_MESH::drawConstrainedNodes()
{
  glDisable(GL_DEPTH_TEST);
    glColor3f(1.0, 0.0, 0.0);
    glPointSize(10.0f);
    glBegin(GL_POINTS);
      for (unsigned int x = _unconstrainedSize; x < _vertices.size(); x++)  
        glVertex3f(_vertices[x][0], _vertices[x][1], _vertices[x][2]);
    glEnd();
  glEnable(GL_DEPTH_TEST);
}

void TET_MESH::drawLowresEmbedding()
{
  // glDisable(GL_DEPTH_TEST);
  //   glColor3f(1.0, 0.0, 0.0);
  //   glPointSize(10.0f);
  //   glBegin(GL_POINTS);
  //     for (unsigned int x = 0; x < _lowresVertices.size(); x++)  
  //       glVertex3f(_lowresVertices[x][0], _lowresVertices[x][1], _lowresVertices[x][2]);
  //   glEnd();
  // glEnable(GL_DEPTH_TEST);
  for (unsigned int x = 0 ; x < _lowresSurfaceFaces.size(); x++){
    _lowresSurfaceFaces[x].draw();
  }
}

void TET_MESH::drawDisplacement()
{
  Real maxDisp = 0;
  for(unsigned int i = 0; i < _x.size() / 3; i++)
  {
    Real disp = _x.segment<3>(i * 3).norm();
    maxDisp = disp > maxDisp ? disp : maxDisp;
  }
  glColor3f(1.0, 1.0, 1.0);
  for (unsigned int x = 0 ; x < _surfaceFaces.size(); x++){
    bool draw = false;
    for(int y = 0; y < 3 && !draw; y++){
      int id = _vertexID[_surfaceFaces[x].vertices[y]];
      Real diff = (_vertices[id] - _restPose[id]).norm();
      // if(diff / maxDisp > 0.03)
      if(diff > 0)
        draw = true;
    }
    if(draw)
      _surfaceFaces[x].draw();
  }
}
void TET_MESH::drawSurfaceVertexOneRing(int surfaceVertexID)
{
  if(surfaceVertexID < 0 || surfaceVertexID >= _surfaceVertices.size())
    return;

  if(_surfaceVertexOneRings.empty())
    computeSurfaceVertexOneRings();

  glDisable(GL_DEPTH_TEST);
    glPointSize(10.0f);
    glBegin(GL_POINTS);
    
      glColor3f(1.0, 0.0, 0.0);
      VEC3F& sv = *_surfaceVertices[surfaceVertexID];
      glVertex3f(sv[0], sv[1], sv[2]);

      glColor3f(0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < _surfaceVertexOneRings[surfaceVertexID].size(); x++){
        VEC3F& node = *_surfaceVertexOneRings[surfaceVertexID][x];
        glVertex3f(node[0], node[1], node[2]);
      }
    glEnd();
  glEnable(GL_DEPTH_TEST);
}
void TET_MESH::drawCollisionPairs()
{
  glColor4f(1.0, 0.0, 0.0, 1.0);
  for(unsigned int i = 0; i < _collisionPairs.size(); i++){
    VEC3F contactPoint = _collisionPairs[i].second->contactPoint(*(_collisionPairs[i].first));
    VEC3F& vertex = *(_collisionPairs[i].first);
    glPointSize(5.0f);
      glBegin(GL_POINTS);
        glColor4f(0.0, 0.0, 1.0, 1.0);
        glVertex3f(vertex[0], vertex[1], vertex[2]);
        glColor4f(1.0, 0.0, 0.0, 1.0);
        glVertex3f(contactPoint[0], contactPoint[1], contactPoint[2]);
      glEnd();

    glLineWidth(2.0f);
      glColor4f(1.0, 1.0, 0.0, 1.0);
      glBegin(GL_LINES);
        glVertex3f(vertex[0], vertex[1], vertex[2]);
        glVertex3f(contactPoint[0], contactPoint[1], contactPoint[2]);
      glEnd();
  }
}
