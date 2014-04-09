/**
 * \file Data.h
 * \brief Data module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to handle data structure
 */

#ifndef DATA_H_
#define DATA_H_

/**
 * \struct Data_structure
 * \brief Data structure
 * 
 * Structure containing all the information about each cluster
 */

// AS ** Check Nfilt.. For the moment is fixed to 30


typedef struct{
  //  char name[20]; /** name of the cluster */
  //  char field[20]; /** SPT field name */
  int clus_id;
  int field_id;  
  double redshift; /** redshift of teh cluster */
  double lambda; /** lambda observable */
  double dlambda; /** lambda observable */
  double y0[30]; /** Y0 */
  double dy0[30]; /** measurement uncertainty on Dy0 */
  double filter[30]; /** Filtering scale of Arnaud profile - R500 arcmin **/
} Data_structure;


typedef struct{
  double r500[146];
  double y0_to_Y500cyl[146];
} y0_to_Y500cyl_structure;

Data_structure* Read_Data(char*, int*);
void Free_Data(Data_structure*);
void Copy_Data(Data_structure, Data_structure*);


#endif
