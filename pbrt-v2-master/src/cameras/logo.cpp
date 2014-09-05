#include "vdb.h"
#include <math.h>

#define PI 3.14159265
#define PARTICLERADIUS 0.4
#define WAVEAMPLITUDE 0.75
#define FREQ 0.75
#define WATERDEPTH 2.0
#define YMAX 10
#define ZMAX 21

void drawParticle(double x, double y, double z, double r, double g, double b){
  int fineness = 200;
  double unit = 2.0 * PI/(double)fineness;
  vdb_color(r,g,b);
  for(int i = 0; i < fineness; i++){
    double start1 = PARTICLERADIUS * cos((double)i * unit); //+ x;
    double start2 = PARTICLERADIUS * sin((double)i * unit); //+ y;
    double end1 = PARTICLERADIUS * cos((double)(i+1) * unit); //+ x; 
    double end2 = PARTICLERADIUS * sin((double)(i+1) * unit); //+ y;
  	vdb_line(x + start1, y + start2, z, x + end1, y + end2, z);
  	vdb_line(x + start1, y, z + start2, x + end1, y, z + end2);
  	vdb_line(x, y + start1, z + start2, x, y + end1, z + end2);
  }	
}

void drawRing(double radius, double x, double z, double r, double g, double b){
	int fineness = 24;
	double unit = 2.0 * PI/(double)fineness;
	for(int i = 0; i < fineness; i++){
		double renderZ = z + radius * cos((double)i * unit);		
		double renderX = x + radius * sin((double)i * unit);
		double renderY = 4.0;
		while(renderY < 6.0){
			drawParticle(renderX, renderY, renderZ, r, g, b);
			renderY += 0.5;
		}
	}
}

void drawWater(double t){
	double palette[12][3] = {  
	 {1,0,0} ,   
	 {1,0.5,0} ,   
	 {1,1,0} ,
	 {0.5,1,0} ,
	 {0,1,0} ,
	 {0,1,0.5} ,
	 {0,1,1} ,
	 {0,0.5,1} ,
	 {0,0,1} ,
	 {0.5,0,1} ,
	 {1,0,1} ,
	 {1,0,0.5}
	};


  double z = 0.0;
  int count = ((int)(t/0.1))%12;
  while(z < ZMAX){
  	vdb_color(palette[count][0], palette[count][1], palette[count][2]);
  	double height = WATERDEPTH + WAVEAMPLITUDE * sin(FREQ*(z + t));
  	double x = 0.0;
  	while(x < height){
  		double y = 0.0;
  		while(y < YMAX){
  			vdb_point(x,y,z);
  			y += 0.5;
  		}
  		x += 0.5;
  	}
  	z += 0.5;
  	count ++;
  	if(count == 12) count = 0;  	
  }
}

int main(){
  double t = 0.0;
  while(true){
  	vdb_frame();
  	vdb_begin();
  	drawWater(t);
  	t += 0.1;
  	double ringRadius = 2.7;
  	double bottom = 5.7;
  	double top = bottom + ringRadius;
  	drawRing(ringRadius, top, 3.6, 1, 1, 0); //top left
  	drawRing(ringRadius, top, 10, 1, 1, 1); //top center
  	drawRing(ringRadius, top, 16.4, 0, 1, 1); //top right
  	drawRing(ringRadius, bottom, 6.8, 0, 0, 1); //bottom left
  	drawRing(ringRadius, bottom, 13.2, 1, 0, 1); //bottom left
  	vdb_end();
  }
  return 0;  
}