#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #define AQUARIUS
#define makeLong(hi, low) (((long) hi) << 32 | (low))
struct Gadget_particle
{

  float Pos[3];
  float Vel[3];
  float Mass;

  int Type;

  float Rho, U, Temp, Ne;

};

const float Mpc2kpc=1000.;
struct Gadget_particle *P;
int *PIDmap;
unsigned int *Id;

struct gadget_io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int nhighword[6];
  int filler[16];
  // char filler[256-6*4-6*8-8-8-4-4-6*4-4-4-8-8-8-8-4-4-6*4];
  // int hashtabsize;
  //char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 11]; /* fills to 256 Bytes */
};

long gadget_load_snapshot(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  long longdummy;
  int  dummy;
  long NumPart;
  int i, j, k, ntot_withmasses,Ngas;
  int t, n, off, pc, pc_new, pc_sph,local_nids;
#ifdef AQUARIUS
  double *tmp;
#else
  float *tmp;
#endif
  struct gadget_io_header header1;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  tmp = malloc(0);
  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      printf("reading `%s' ...\n", buf);
      fflush(stdout);

      // fread(&dummy, sizeof(dummy), 1, fd);
      SKIP;
      printf("dummy = %d\n",dummy);
      fread(&header1, sizeof(header1), 1, fd);
      fseek(fd,256+sizeof(dummy),SEEK_SET);
      SKIP;
      // fread(&dummy, sizeof(dummy), 1, fd);
      // printf("dummy = %d\n",dummy);
      //SKIP;
      if(files == 1)
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    {
	      NumPart += header1.npart[k];
	    }
	  Ngas = header1.npart[0];
	}
      else
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    {
	      NumPart += makeLong(header1.nhighword[k],header1.npartTotal[k]);
	      printf("Mass[%d] = %lf\n",k,header1.mass[k]);
	      printf("Total num %d : %d\n",k,header1.npartTotal[k]);
	      printf("high word %d : %d\n",k,header1.nhighword[k]);
	    }
	  Ngas = header1.npartTotal[0];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	{
	  printf("Total numpart : %ld\n",NumPart);
	  
	  //allocate_memory();
	  /* P = malloc(NumPart * sizeof(struct Gadget_particle)); */
	  /* Id = malloc(NumPart * sizeof(int)); */
	  /* PIDmap = malloc(NumPart * sizeof(int)); */
	  /* P--; */
	  /* Id--; */
	  /* PIDmap--; */
	}
      local_nids = 0;
      for(k=0;k<6;k++)
	{
	  printf("N[%d] : %d\n",k,header1.npart[k]);
	  local_nids += header1.npart[k];
	}
      /* for(k=0;k<6;k++) */
      /* 	{ */
      /* 	  printf("Total N[%d] : %d\n",k,header1.npartTotal[k]); */
      /* 	} */
      /* for(k=0;k<6;k++) */
      /* 	{ */
      /* 	  printf("mass[%d] : %lf\n",k,(double)header1.mass[k]); */
      /* 	} */

      printf("time : %lf\n",header1.time);
      printf("redshift : %lf\n",header1.redshift);
      printf("flag_sfr : %d\n", header1.flag_sfr);
      printf("flag_feedback : %d\n", header1.flag_feedback);
      printf("flag_cooling : %d\n", header1.flag_cooling);
      printf("numfiles : %d\n", header1.num_files);
      printf("BoxSize : %lf\n",header1.BoxSize);
      printf("Omega0 : %lf\n",header1.Omega0);
      printf("OmegaLambda : %lf\n",header1.OmegaLambda);
      printf("HubbleParam : %lf\n",header1.HubbleParam);
      for(k=0;k<16;k++)
      	{
      	  printf("fillers[%d] : %d\n",k, header1.filler[k]);
      	}
      
      SKIP;
      // fread(&dummy, sizeof(dummy), 1, fd);
      printf("dummy = %d/%d\n",dummy, local_nids);
#ifdef AQUARIUS
      tmp = realloc(tmp,sizeof(double)*3);
#else
      tmp = realloc(tmp,sizeof(float)*3);
#endif
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
#ifdef AQUARIUS
	      fread(&tmp[0], sizeof(double), 3, fd);
#else
	      fread(&tmp[0], sizeof(float), 3, fd);
#endif
	      if(n%100000 == 0) printf("%lf %lf %lf\n",tmp[0],tmp[1],tmp[2]);
	      /* for(j=0;j<3;j++) */
	      /* 	{ */
	      /* 	  printf("%lf %lf %lf\n",tmp[0],tmp[1],tmp[2]); */
	      /* 	} */

	      pc_new++;
	    }
	}
      SKIP; 
      // fread(&dummy, sizeof(dummy), 1, fd);
      printf("dummy = %d\n",dummy);

      SKIP;

      // fread(&dummy, sizeof(dummy), 1, fd);
      printf("dummy = %d/%d\n",dummy, local_nids);
#ifdef AQUARIUS
      tmp = realloc(tmp,sizeof(double)*3);
#else
      tmp = realloc(tmp,sizeof(float)*3);
#endif
     for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
#ifdef AQUARIUS
	      fread(&tmp[0], sizeof(double), 3, fd);
#else
	      fread(&tmp[0], sizeof(float), 3, fd);
#endif

	      //for(j=0;j<3;j++)
	        //printf("%f %f %f\n",tmp[0],tmp[1],tmp[2]);
		// P[pc_new].Vel[j] = (float) tmp[j]*(header1.redshift+1.);
	      pc_new++;
	    }
	}
      SKIP; 
      // fread(&dummy, sizeof(dummy), 1, fd);
      printf("dummy = %d\n",dummy);

      SKIP;
      // fread(&dummy, sizeof(dummy), 1, fd);
      printf("dummy = %d/%d\n",dummy, local_nids);
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      // fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP; 
      // fread(&dummy, sizeof(dummy), 1, fd);
      printf("dummy = %d\n",dummy);


      if(ntot_withmasses > 0)
	{
	  SKIP;
	  // fread(&dummy, sizeof(dummy), 1, fd);
	  printf("dummy = %d/%d\n",dummy, local_nids);
#ifdef AQUARIUS
	  tmp = realloc(tmp,sizeof(double));
#else
	  tmp = realloc(tmp,sizeof(float));
#endif
	  /* fread(&dummy, sizeof(dummy), 1, fd); */
	  /* printf("dummy = %d\n",dummy); */
	}
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      //  P[pc_new].Type = k;

	      if(header1.mass[k] == 0)
		{
#ifdef AQUARIUS
		  fread(&tmp[0], sizeof(double), 1, fd);
#else
		  fread(&tmp[0], sizeof(float), 1, fd);
#endif
		  // P[pc_new].Mass = (float) tmp[0];
		}
	      else
		// P[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	{
	  SKIP;
	  // fread(&dummy, sizeof(dummy), 1, fd);
	  printf("dummy = %d\n",dummy);
	}

      if(header1.npart[0] > 0)
	{
	  SKIP;
	  // fread(&dummy, sizeof(dummy), 1, fd);
	  printf("dummy = %d\n",dummy);
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      // fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP; 
	  // fread(&dummy, sizeof(dummy), 1, fd);
	  printf("dummy = %d\n",dummy);
	  SKIP;
	  // fread(&dummy, sizeof(dummy), 1, fd);
	  printf("dummy = %d\n",dummy);
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      // fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP; 
	  // fread(&dummy, sizeof(dummy), 1, fd);
	  printf("dummy = %d\n",dummy);	  
	  if(header1.flag_cooling)
	    {
	      // SKIP; 
	      // fread(&dummy, sizeof(dummy), 1, fd);
	      printf("dummy = %d\n",dummy);
	      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
		{
		  // fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      // SKIP;
	      // fread(&dummy, sizeof(dummy), 1, fd);
	      printf("dummy = %d\n",dummy);
	      /* SKIP; */
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	      {
		// P[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }
  //printf("%d %f %f %f\n",Id[NumPart],P[NumPart].Pos[0],P[NumPart].Vel[0],P[NumPart].Mass);
  /* for(i=1;i<=NumPart;i++) */
  /*   { */
  /*     //if(P[i].Type != 1) printf("%d %d %f\n",(int)i,(int)Id[i],P[i].Mass); */
  /*     //printf("%d => %d\n",i,(int)Id[i]); */
  /*     //if(i != Id[i]) printf("%d => %d\n",i,(int)Id[i]); */
  /*     PIDmap[Id[i]] = i; */
  /*   } */
  //printf("%d %f %f %f\n",Id[NumPart],P[NumPart].Pos[0],P[NumPart].Vel[0],P[NumPart].Mass);
  return NumPart;
}

int main ()
{
  printf("%ld\n",gadget_load_snapshot("/scratch/01937/cs390/cubepm_130315_6_1728_47Mpc_ext2/snapdir_000/cube2gadget_000",216));
}
