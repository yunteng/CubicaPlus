/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
template<class T>
void VIEWER<T>::init()
{
	Real moveSpeed = 1.0;
  Real eyeDistanceScale = 1.0;

  // moveSpeed = configFile.getFloat("navigation speed", moveSpeed);
  // eyeDistanceScale = configFile.getFloat("view distance", eyeDistanceScale);

  glvu.SetMoveSpeed( moveSpeed * glvu.GetMoveSpeed() );

  int argc = 0;
  glutInit(&argc, NULL);

  glvuVec3f boxCenter(0.5, 0.5, 0.5);
  glvuWindow( boxCenter, eyeDistanceScale );
}
template<class T>
void VIEWER<T>::displayFunc()
{
	static GLfloat mat_shininess[] = { 80.0 };

	// static GLfloat mat_diffuse[] = { 0.5, 0.2, 0.3, 1.0 };
  static GLfloat mat_diffuse[] = { 0.96, 0.96, 0.86, 1.0 };
	static GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };

  Camera* camera = glvu.GetCurrentCam();
  glvuVec3f Eye, Lookat, Up;
  camera->GetLookAtParams(&Eye, &Lookat, &Up);

  // repack camera settings into arrays
  float eye[] = {Eye.x, Eye.y, Eye.z};
  float look[3];
  look[0] = Lookat.x - eye[0];
  look[1] = Lookat.y - eye[1];
  look[2] = Lookat.z - eye[2];
  float magnitude = 1.0f / sqrt(look[0] * look[0] + look[1] * look[1] + look[2] * look[2]);
  look[0] *= magnitude;
  look[1] *= magnitude;
  look[2] *= magnitude;

  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    
    // glColor4f(0.05, 0.75, 1.0, 1.0);
    glColor4f(0.96, 0.96, 0.86, 1.0);
    // glScalef(1, 1, 1);
    // glutWireCube(1.0);
    simulator->display();

  glvu.EndFrame();
}

template<class T>
void VIEWER<T>::keyboardFunc(unsigned char Key, int x, int y)
{
  simulator->keyboardFunc(Key);
  switch(Key)
  {
  	case 'a':{
  		animate = !animate;
  		break;
  	}
    case 's':{
      step = !step;
      break;
    }
    // case 'Z':{
      
    //   break;
    // }
    case 'v':
      {
        Camera* camera = glvu.GetCurrentCam();
        glvuVec3f eye;
        glvuVec3f lookat;
        glvuVec3f up;
        camera->GetLookAtParams(&eye, &lookat, &up);
        cout << " Eye(" << eye[0] << ", " << eye[1] << ", " << eye[2] << "), " ;
        cout << " LookAtCntr(" << lookat[0] << ", " << lookat[1] << ", " << lookat[2] << "), " ;
        cout << " Up(" << up[0] << ", " << up[1] << ", " << up[2] << ");\n";

        cout << "eye = " << eye[0] << "," << eye[1] << "," << eye[2]
        		 << "\nlookatcntr = " << lookat[0] << "," << lookat[1] << "," << lookat[2]
        		 << "\nup = " << up[0] << "," << up[1] << "," << up[2] << "\n";

      }
      break;
    case 'q':
    case 'Q':
      TIMING_BREAKDOWN::printTimingBreakdown();
      // TIMING_BREAKDOWN::writeFrameTime("./frametime.txt");
      exit(0);
      break;
  };

  glutPostRedisplay();
  if (Key != '=')
    glvu.Keyboard(Key,x,y);
}

template<class T>
void VIEWER<T>::idleFunc()
{
	if(animate){
    bool getScreenshot = SIMPLE_PARSER::getBool("get screen shot", false);
    if(getScreenshot)
      screenshot(simulator->renderPath, simulator->previousFrame);
		simulator->step();
    if(step)
      animate = false;
    glutPostRedisplay(); 
  }
}

template<class T>
void VIEWER<T>::mouseFunc(int button, int state, int x, int y)
{
  static float clickZ;
  int Modifiers = glutGetModifiers();
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_DOWN &&
      Modifiers & GLUT_ACTIVE_SHIFT)
  {
    // retrieve and store the depth of this click
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glReadPixels(x, viewport[3] - y, 
                 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &clickZ);

    // get the world space coordinate
    VEC3F point = unproject(x,y,clickZ);

    // hand the coordinate to the tet mesh
    // simulator->click(point);
    // mouseClicked = true;
    return;
  }
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_UP)
  {
    // mouseClicked = false;
    // simulator->unclick();
    // integrator->unclick();
    return;
  } 
	// pass through to default handler
  glvu.Mouse(button,state,x,y);
}

template<class T>
void VIEWER<T>::motionFunc(int x, int y)
{
	glvu.Motion(x,y);
}

template<class T>
int VIEWER<T>::glvuWindow(glvuVec3f bboxCenter, Real eyeDistanceScale)
{
	glvu.Init((char *)"GLVU Window",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            windowStartX, windowStartY, windowWidth, windowHeight);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.8, 0.8, 0.8, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  //glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glShadeModel(GL_SMOOTH);
  //glShadeModel(GL_FLAT);
  // glClearColor(1, 1, 1, 1);
  glClearColor(0, 0, 0, 1);

  glutDisplayFunc(VIEWER<T>::displayFunc);
  glutMouseFunc(VIEWER<T>::mouseFunc);
  glutMotionFunc(VIEWER<T>::motionFunc);
  glutKeyboardFunc(VIEWER<T>::keyboardFunc);
  glutIdleFunc(VIEWER<T>::idleFunc);

  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10);

  VEC3F myEye = SIMPLE_PARSER::getVEC3F("eye", VEC3F(0, 0.3, 3));
  glvuVec3f Eye(myEye[0], myEye[1], myEye[2]);

  VEC3F myLookAtCntr = SIMPLE_PARSER::getVEC3F("lookatcntr", VEC3F(0, 0.2, 2));
  glvuVec3f LookAtCntr(myLookAtCntr[0], myLookAtCntr[1], myLookAtCntr[2]);

  VEC3F myUp = SIMPLE_PARSER::getVEC3F("up", VEC3F(0, 0.9, -0.2));
  glvuVec3f Up(myUp[0], myUp[1], myUp[2]);
  
  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov,Aspect, Near, Far);

  // reset center
  glvuVec3f center(0.5f, 0.5f, 0.5f);
  glvu.SetWorldCenter(center);
  
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
// Dump a frame with this number to this directory
//////////////////////////////////////////////////////////////////////////////
template<class T>
void VIEWER<T>::screenshot(string renderPath, int frame)
{
  FILE *fp;

  string number = IO::itoPaddedString(frame);
  // char FileName[256];
  string filename = renderPath + SIMPLE_PARSER::getString("render prefix", "renderGL") + "." + number + ".ppm";
  // sprintf(FileName,"%srenderGL.%s.ppm", renderPath.c_str(), number.c_str());

  GLint OldReadBuffer;
  glGetIntegerv(GL_READ_BUFFER,&OldReadBuffer);
  glReadBuffer(GL_FRONT);

  GLint OldPackAlignment;
  glGetIntegerv(GL_PACK_ALIGNMENT,&OldPackAlignment); 
  glPixelStorei(GL_PACK_ALIGNMENT,1);

  int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
  int NumPixels = WW*WH;
  GLubyte* Pixels = new GLubyte[NumPixels*3];
  if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
  glReadPixels(0,0,WW,WH,GL_RGB,GL_UNSIGNED_BYTE,Pixels);

  fp = fopen(filename.c_str(), "wb");
  fprintf(fp, "P6\n%d %d\n255\n", WW, WH);
  for(int y = WH - 1; y >= 0; y--){
    fwrite(Pixels + y * WW * 3, 1, WW * 3, fp);
  }
  // fwrite(Pixels,1,NumPixels*3,fp);
  fclose(fp);
  delete[] Pixels;

  glPixelStorei(GL_PACK_ALIGNMENT,OldPackAlignment);
  glReadBuffer((GLenum)OldReadBuffer);
}

//////////////////////////////////////////////////////////////////////////////
// Translate screen space to world space
//
// Adapted from the NeHe page:
// http://nehe.gamedev.net/data/articles/article.asp?article=13
//////////////////////////////////////////////////////////////////////////////
template<class T>
VEC3F VIEWER<T>::unproject(float x, float y, float z)
{
  GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

  double worldX, worldY, worldZ;
	gluUnProject(x, viewport[3] - y, z,
               modelview, projection, viewport, 
               &worldX, &worldY, &worldZ);

  return VEC3F(worldX, worldY, worldZ);
}


