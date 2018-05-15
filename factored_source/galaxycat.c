#include "galaxycat.h"

GalaxyCatalog *read_galaxies_from_file(char *filename, int npart);
int  count_lines_in_file(char *filename);
void outputGalaxies(GalaxyCatalog *cat, char *filename);
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box);
void free_cat(GalaxyCatalog *cat);

//====================================================
// Read galaxies from file. Format: RA DEC z
//====================================================
GalaxyCatalog *read_galaxies_from_file(char *filename, int ngalaxies){
  int i;

  printf("\n====================================\n"); 
  printf("Reading from file [%s]\n", filename); 
  printf("====================================\n"); 
  printf("Galaxy file has %d galaxies\n",ngalaxies);

  // Allocate particle array
  GalaxyCatalog *cat = malloc(sizeof(GalaxyCatalog));
  cat->galaxies = malloc(sizeof(Galaxy)*ngalaxies);
  Galaxy *allgalaxies = cat->galaxies;
  cat->ngalaxies = ngalaxies;
  cat->allocated = 1;
    
  double sum_w = 0.0, sum_w2 = 0.0;

  // Read the data
  FILE *fp = fopen(filename, "r");
  for(i = 0; i < ngalaxies; i++){
    char line[1024];
    double Pos[3], RA, DEC, z, w = 1.0;

    // Read galaxy
    fgets(line, sizeof(line),fp);
#ifdef WEIGHTS
    sscanf(line, "%lf %lf %lf %lf", &RA, &DEC, &z, &w);
#else
    sscanf(line, "%lf %lf %lf", &RA, &DEC, &z);
#endif

    // Convert to positions in Mpc/h
    double costheta = cos(2.0 * M_PI / 360.0 * (90.0 - DEC));
    double sintheta = sqrt(1.0 - costheta*costheta);
    double phi      = RA * 2.0 * M_PI / 360.0;
    double r        = r_of_z(z);

    // Set positions
    Pos[0] = r * sintheta * cos(phi);
    Pos[1] = r * sintheta * sin(phi);
    Pos[2] = r * costheta;

    // Store galaxies
    allgalaxies[i].x[0] = Pos[0];
    allgalaxies[i].x[1] = Pos[1];
    allgalaxies[i].x[2] = Pos[2];
#ifdef WEIGHTS
    allgalaxies[i].w = w;
#endif  
    sum_w  += w;
    sum_w2 += w*w;
  }
  fclose(fp);
  
  // The mean weight and RMS
  printf("Mean weight: %lf  RMS: %lf\n", sum_w/(double)ngalaxies, sqrt(sum_w2/(double)ngalaxies));
  cat->sum_w = sum_w;
  cat->sum_w2 = sum_w2;

  return cat;
}

//====================================================
// Compute the boxsize we need to encompas all galaxies
// in both catalogs and shift the positions so that 
// they are inside [0, BOX]^3
//====================================================
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box){
  int ngalaxies  = cat->ngalaxies;
  Galaxy *galaxies = cat->galaxies;
  int ngalaxies2 = cat2->ngalaxies;
  Galaxy *galaxies2 = cat2->galaxies;
  
  // Compute max and min position
  double max_x = -1e100, min_x = 1e100;
  double max_y = -1e100, min_y = 1e100;
  double max_z = -1e100, min_z = 1e100;
 
  int i;
  for(i = 0; i < ngalaxies; i++){
    double *Pos = &galaxies[i].x[0];
    if(Pos[0] > max_x) max_x = Pos[0];
    if(Pos[1] > max_y) max_y = Pos[1];
    if(Pos[2] > max_z) max_z = Pos[2];
    if(Pos[0] < min_x) min_x = Pos[0];
    if(Pos[1] < min_y) min_y = Pos[1];
    if(Pos[2] < min_z) min_z = Pos[2];
  }
  for(i = 0; i < ngalaxies2; i++){
    double *Pos = &galaxies2[i].x[0];
    if(Pos[0] > max_x) max_x = Pos[0];
    if(Pos[1] > max_y) max_y = Pos[1];
    if(Pos[2] > max_z) max_z = Pos[2];
    if(Pos[0] < min_x) min_x = Pos[0];
    if(Pos[1] < min_y) min_y = Pos[1];
    if(Pos[2] < min_z) min_z = Pos[2];
  }

  // Shift positions
  for(i = 0; i < ngalaxies; i++){
    double *Pos = &galaxies[i].x[0];
    Pos[0] -= min_x;
    Pos[1] -= min_y;
    Pos[2] -= min_z;
  }
  for(i = 0; i < ngalaxies2; i++){
    double *Pos = &galaxies2[i].x[0];
    Pos[0] -= min_x;
    Pos[1] -= min_y;
    Pos[2] -= min_z;
  }

  // New max_x,y,z values
  max_x = max_x-min_x;
  max_y = max_y-min_y;
  max_z = max_z-min_z;

  // The min/max positions (separations)
  printf("\n====================================\n");
  printf("Shifting particles and computing boxsize:\n");
  printf("====================================\n");
  printf("Min/Max X position: 0.0 -> %5.2lf Mpc/h\n", max_x);
  printf("Min/Max Y position: 0.0 -> %5.2lf Mpc/h\n", max_y);
  printf("Min/Max Z position: 0.0 -> %5.2lf Mpc/h\n", max_z);
  
  // The largest displacement in either direction
  *box = 1.01 * MAX(max_x, MAX(max_y, max_z));
  printf("The boxsize we will use is %5.2lf Mpc/h\n", *box);
}

//====================================================
// Count the number of lines in a file
//====================================================
int count_lines_in_file(char *filename){
  int n = 0;
  FILE *fp = fopen(filename, "r");
  while(!feof(fp)){
    char ch = fgetc(fp);
    if(ch == '\n') n++;
  }
  return n;
}

//====================================================
// Free the memory associated with a galaxy catalog
//====================================================
void free_cat(GalaxyCatalog *cat){
  if(cat->allocated){
    free(cat->galaxies);
    free(cat);
  }
}

//====================================================
// Output galaxy catalog in physical coordinates
//====================================================
void outputGalaxies(GalaxyCatalog *cat, char *filename){
  int ngalaxies = cat->ngalaxies;
  Galaxy *galaxies = cat->galaxies;
  FILE *fp = fopen(filename, "w");

  int i;
  for(i = 0; i < ngalaxies; i++){
    Galaxy *curgalaxy = &galaxies[i];
#ifdef WEIGHTS
    fprintf(fp, "%lf  %lf  %lf  %lf\n", curgalaxy->x[0], curgalaxy->x[1], curgalaxy->x[2], curgalaxy->w);
#else
    fprintf(fp, "%lf  %lf  %lf\n", curgalaxy->x[0], curgalaxy->x[1], curgalaxy->x[2]);
#endif
  }
  fclose(fp);
}

