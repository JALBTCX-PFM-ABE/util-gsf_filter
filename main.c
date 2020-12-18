
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include "gsf_filter.h"
#include "version.h"


extern int32_t gsfError;


void usage ()
{
      fprintf (stderr, "USAGE: gsf_filter [--std STD] [--deep] GSF_FILE\n\n");
      fprintf (stderr, "Where:\n");
      fprintf (stderr, "\tGSF_FILE = Path to GSF file.\n");
      fprintf (stderr, "\tSTD = Optional number of standard deviations to filter (default = 2.0)\n");
      fprintf (stderr, "\t-d = Filter only in the downward (deep filter) direction\n\n");
}


int32_t main (int32_t argc, char **argv)
{
  gsfDataID           id;
  gsfRecords          gsf_record;
  int32_t             hnd, i, j, k, recnum, percent = 0, old_percent = -1, ret, start_rec, count, page_size = 1000;
  int32_t             grid_height = 0, grid_width = 0, xn, yn, *aping = NULL, prev_ping = -1, option_index = 0;
  int16_t             *abeam = NULL;
  float               std_env, dep, *adep = NULL, avg_z;
  double              *alat = NULL, *alon = NULL, lateral, ang1, ang2, lat, lon, sum_z, sum2_z, grid_size;
  double              dx, rlat1, rlat2, rlon1, rlon2, az, prev_lat = -999.0, prev_lon = -999.0;
  NV_F64_COORD2       xy2, nxy;
  NV_F64_XYMBR        mbr;
  char                c, comment[16384], file[512];
  uint8_t             endloop = NVFalse, skipflag = NVFalse, flushflag = NVFalse, deepflag = NVFalse;
  GRID_REC            **grid = NULL;
  extern char         *optarg;
  extern int          optind;
  static struct option long_options[] = {{"std", required_argument, 0, 0},
                                         {"deep", no_argument, 0, 0},
                                         {0, no_argument, 0, 0}};



  int32_t write_history (int32_t, char **, char *, char *, int32_t);


  printf ("\n\n %s \n\n",VERSION);


  deepflag = NVFalse;
  std_env = 2.0;


  while (NVTrue) 
    {
      c = (char) getopt_long (argc, argv, "d", long_options, &option_index);
      if (c == -1) break;

      switch (c) 
        {
        case 0:

          switch (option_index)
            {
            case 0:
              sscanf (optarg, "%f", &std_env);
              if (std_env < 1.0 || std_env > 10.0) std_env = 2.0;
              break;

            case 1:
              deepflag = NVTrue;
              break;
            }
          break;

        default:
          usage ();
          exit (-1);
          break;
        }
    }


  /* Make sure we got the mandatory file name argument.  */

  if (optind >= argc)
    {
      usage ();
      exit (-1);
    }


  strcpy (file, argv[optind]);

  if (gsfOpen (file, GSF_UPDATE_INDEX, &hnd))
    {
      gsfPrintError (stderr);
      exit (-1);
    }


  printf ("File : %s\n\n", file);


  start_rec = recnum = 1;


  while (!endloop)
    {
      count = 0;
      mbr.min_x = 999.0;
      mbr.max_x = -999.0;
      mbr.min_y = 999.0;
      mbr.max_y = -999.0;
      sum_z = 0.0;
      skipflag = NVFalse;


      /*  Read "page_size" pings and load them into local memory.  */

      for (j = start_rec ; j < start_rec + page_size ; j++)
        {
          id.recordID = GSF_RECORD_SWATH_BATHYMETRY_PING;
          id.record_number = j;

          if (gsfRead (hnd, GSF_RECORD_SWATH_BATHYMETRY_PING, &id, &gsf_record, NULL, 0) < 0)
            {
              endloop = NVTrue;

              if (gsfError != GSF_INVALID_RECORD_NUMBER) gsfPrintError(stderr);
              break;
            }


          lat = gsf_record.mb_ping.latitude;
          lon = gsf_record.mb_ping.longitude;


          /*  Only deal with valid pings.  */

          if ((lat <= 90.0) && (lon <= 180.0) && !(gsf_record.mb_ping.ping_flags & GSF_IGNORE_PING))
            {
              /*  If we jumped more than 1000 meters we want to close this box, process it,
                  and then start over with this record.  */

              invgp (NV_A0, NV_B0, lat, lon, prev_lat, prev_lon, &dx, &az);
              if (prev_lat > -900.0 && dx > 1000.0)
                {
                  start_rec = j;
                  skipflag = NVTrue;
                  break;
                }
              prev_lat = lat;
              prev_lon = lon;


              ang1 = gsf_record.mb_ping.heading + 90.0;
              ang2 = gsf_record.mb_ping.heading;


              for (i = 0 ; i < gsf_record.mb_ping.number_beams ; i++) 
                {
                  dep = gsf_record.mb_ping.depth[i];
                  if (dep == 0.0 && gsf_record.mb_ping.nominal_depth != NULL) dep = gsf_record.mb_ping.nominal_depth[i];


                  /*  Only deal with valid beams.  */

                  if (dep != 0.0 && gsf_record.mb_ping.beam_flags != NULL && 
                      !(check_flag (gsf_record.mb_ping.beam_flags[i], NV_GSF_IGNORE_NULL_BEAM)) &&
                      !(check_flag (gsf_record.mb_ping.beam_flags[i], (NV_GSF_IGNORE_MANUALLY_EDITED | 
                                                                       NV_GSF_IGNORE_FILTER_EDITED)))) 
                    {
                      /*  Adjust for cross track position.  */

                      lateral = gsf_record.mb_ping.across_track[i];
                      newgp (lat, lon, ang1, lateral, &nxy.y, &nxy.x);


                      /*  if the along track array is present then use it  */
                
                      if (gsf_record.mb_ping.along_track != (double *) NULL) 
                        {
                          xy2.y = nxy.y;
                          xy2.x = nxy.x;
                          lateral = gsf_record.mb_ping.along_track[i];

                          newgp (xy2.y, xy2.x, ang2, lateral, &nxy.y, &nxy.x);
                        }


                      /*  Reallocate the point memory.  */

                      alat = (double *) realloc (alat, (count + 1) * sizeof (double));
                      if (alat == NULL)
                        {
                          perror ("Allocating alat memory");
                          exit (-1);
                        }

                      alon = (double *) realloc (alon, (count + 1) * sizeof (double));
                      if (alon == NULL)
                        {
                          perror ("Allocating alon memory");
                          exit (-1);
                        }

                      adep = (float *) realloc (adep, (count + 1) * sizeof (float));
                      if (adep == NULL)
                        {
                          perror ("Allocating adep memory");
                          exit (-1);
                        }

                      aping = (int32_t *) realloc (aping, (count + 1) * sizeof (int32_t));
                      if (aping == NULL)
                        {
                          perror ("Allocating aping memory");
                          exit (-1);
                        }

                      abeam = (int16_t *) realloc (abeam, (count + 1) * sizeof (int16_t));
                      if (abeam == NULL)
                        {
                          perror ("Allocating abeam memory");
                          exit (-1);
                        }


                      /*  Save the point and compute the mbr.  */

                      alat[count] = nxy.y;
                      alon[count] = nxy.x;
                      adep[count] = dep;
                      aping[count] = j;
                      abeam[count] = i;

                      if (alat[count] < mbr.min_y) mbr.min_y = alat[count];
                      if (alat[count] > mbr.max_y) mbr.max_y = alat[count];
                      if (alon[count] < mbr.min_x) mbr.min_x = alon[count];
                      if (alon[count] > mbr.max_x) mbr.max_x = alon[count];
                      sum_z += adep[count];

                      count++;
                    }
                }
            }
        }


      /*  If we got some points, process them.  */

      if (count)
        {
          /*  Average depth.  */

          avg_z = (float) sum_z / (float) count;


          /*  Compute the grid size.  This is based on a straight down, one-degree footprint size for the
              average depth, times 4.  */

          grid_size = avg_z * 0.017453736 * 4.0 / 111120.0;

          grid_height = NINT (((mbr.max_y - mbr.min_y)) / grid_size + 1.0);
          grid_width = NINT (((mbr.max_x - mbr.min_x)) / grid_size + 1.0);


          /*  Make sure we don't have too large a grid.  */

          while (grid_height > 5000 || grid_width > 5000)
            { 
              grid_size *= 2.0;
              grid_height = NINT (((mbr.max_y - mbr.min_y)) / grid_size + 1.0);
              grid_width = NINT (((mbr.max_x - mbr.min_x)) / grid_size + 1.0);
            }


          /*  Compute the diagonal in meters of a grid cell at the center of the grid.  */

          rlat1 = mbr.min_y + (mbr.max_y - mbr.min_y) / 2.0;
          rlon1 = mbr.min_x + (mbr.max_x - mbr.min_x) / 2.0;
          rlat2 = rlat1 + grid_size;
          rlon2 = rlon1 + grid_size;

          invgp (NV_A0, NV_B0, rlat1, rlon1, rlat2, rlon2, &dx, &az);


          /*  Allocate the grid memory.  */

          grid = (GRID_REC **) calloc (grid_height, sizeof (GRID_REC *));
          if (grid == NULL)
            {
              perror ("Allocating grid memory");
              exit (-1);
            }


          for (i = 0 ; i < grid_height ; i++)
            {
              grid[i] = (GRID_REC *) calloc (grid_width, sizeof (GRID_REC));
              if (grid[i] == NULL)
                {
                  perror ("Allocating grid element memory");
                  exit (-1);
                }


              for (j = 0 ; j < grid_width ; j++) 
                {
                  grid[i][j].count = 0;
                  grid[i][j].depths = NULL;
                }
            }


          /*  Load the grid data from the input points.  */

          for (i = 0 ; i < count ; i++)
            {
              xn = (int32_t) ((alon[i] - mbr.min_x) / grid_size);
              yn = (int32_t) ((alat[i] - mbr.min_y) / grid_size);


              grid[yn][xn].depths = (DEPTH_REC *) realloc (grid[yn][xn].depths, (grid[yn][xn].count + 1) * sizeof (DEPTH_REC));
              if (grid[yn][xn].depths == NULL)
                {
                  perror ("Allocating grid depth memory");
                  exit (-1);
                }

              grid[yn][xn].depths[grid[yn][xn].count].index = i;
              grid[yn][xn].depths[grid[yn][xn].count].filtered = NVFalse;
              grid[yn][xn].count++;
            }


          /*  Free the local position data.  */

          free (alat);
          alat = NULL;
          free (alon);
          alon = NULL;


          /*  Compute the average and standard deviation for each grid node that has data.  */

          for (i = 0 ; i < grid_height ; i++)
            {
              for (j = 0 ; j < grid_width ; j++)
                {
                  if (grid[i][j].count)
                    {
                      grid[i][j].cleared = NVFalse;
                      sum_z = 0.0;
                      sum2_z = 0.0;
                      for (k = 0 ; k < grid[i][j].count ; k++)
                        {
                          sum_z += adep[grid[i][j].depths[k].index];
                          sum2_z += (adep[grid[i][j].depths[k].index] * adep[grid[i][j].depths[k].index]);
                        }

                      grid[i][j].avg = sum_z / (double) grid[i][j].count;

                      if (grid[i][j].count > 1)
                        {
                          grid[i][j].std = sqrt ((sum2_z - ((double) grid[i][j].count * 
                                                            (pow ((double) grid[i][j].avg, 2.0)))) / 
                                                 ((double) grid[i][j].count - 1.0));
                        } 
                      else 
                        {
                          grid[i][j].std = 0.0;
                        }
                    }
                }
            }


          /*  Filter the grid.  */

          gsf_filter (grid, grid_height, grid_width, adep, dx, std_env, deepflag);


          percent = gsfPercent (hnd);
          if (old_percent != percent)
            {
              printf ("%3d%% processed    \r", percent);
              fflush (stdout);
              old_percent = percent;
            }


          /*  Transfer the filtered flags to the adep array by setting the depth to -999999.0 if it is filtered.  
              This way we can do the writes in sequential order instead of bouncing all over the GSF file.  */

          for (i = 0 ; i < grid_height ; i++)
            {
              for (j = 0 ; j < grid_width ; j++)
                {
                  for (k = 0 ; k < grid[i][j].count ; k++)
                    {
                      if (grid[i][j].depths[k].filtered)
                        {
                          adep[grid[i][j].depths[k].index] = -999999.0;
                        }
                    }
                }
            }


          /*  Loop through the data and write out the changed records to the GSF file.  */

          prev_ping = -1;
          for (i = 0 ; i < count ; i++)
            {
              if (adep[i] == -999999.0)
                {
                  if (aping[i] != prev_ping)
                    {
                      if (prev_ping != -1)
                        {
                          id.recordID = GSF_RECORD_SWATH_BATHYMETRY_PING;
                          id.record_number = prev_ping;

                          if (gsfWrite(hnd, &id, &gsf_record) < 0)
                            {
                              gsfPrintError(stderr);
                              exit(-1);
                            }
                          flushflag = NVFalse;
                        }

                      id.recordID = GSF_RECORD_SWATH_BATHYMETRY_PING;
                      id.record_number = aping[i];

                      if (gsfRead (hnd, GSF_RECORD_SWATH_BATHYMETRY_PING, &id, &gsf_record, NULL, 0) < 0)
                        {
                          gsfPrintError (stderr);
                          exit (-1);
                        }
                      prev_ping = aping[i];
                    }
 
                  gsf_record.mb_ping.beam_flags[abeam[i]] |= NV_GSF_IGNORE_FILTER_EDITED;
                  flushflag = NVTrue;
                }
            }
        }


      /*  If we need to flush the last ping, write it out.  */

      if (flushflag)
        {
          id.recordID = GSF_RECORD_SWATH_BATHYMETRY_PING;
          id.record_number = prev_ping;
          if (gsfWrite(hnd, &id, &gsf_record) < 0)
            {
              gsfPrintError (stderr);
              exit (-1);
            }
        }


      /*  Free the depth, ping, and beam memory.  */

      free (adep);
      adep = NULL;
      free (aping);
      aping = NULL;
      free (abeam);
      abeam = NULL;


      /*  Free the grid memory.  */

      if (count)
        {
          for (i = 0 ; i < grid_height ; i++)
            {
              for (j = 0 ; j < grid_width ; j++)
                {
                  if (grid[i][j].count) free (grid[i][j].depths);
                }
              free (grid[i]);
            }
          free (grid);
        }
            

      if (!skipflag) start_rec += page_size;
    }


  percent = gsfPercent (hnd);
  printf ("%3d%% processed    \n", percent);
  gsfClose(hnd);
  printf("\n");
         


  /*  Open the file non-indexed so that we can write a history record.  */

  if (gsfOpen (file, GSF_UPDATE, &hnd))
    {
      gsfPrintError (stderr);
      exit (-1);
    }


  /*  Write a history record describing the filter process.  */

  if (deepflag)
    {
      sprintf (comment, 
               "This file was statistically filtered using the following program and arguments:\n%s -std %.1f -d %s\n",
               argv[0], std_env, file);
    }
  else
    {
      sprintf (comment, 
               "This file was statistically filtered using the following program and arguments:\n%s -std %.1f %s\n",
               argv[0], std_env, file);
    }

  ret = write_history (argc, argv, comment, file, hnd);
  if (ret)
    {
      fprintf(stderr, "Error: %d - writing gsf history record\n", ret);
    }

  gsfClose (hnd);


  return (0);
}
