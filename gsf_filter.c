
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

void gsf_filter (GRID_REC **grid, int32_t height, int32_t width, float *adep, double dx,
                 float std_env, uint8_t deep)
{
  int32_t i, j, m, n, filtered_count, sumcount;
  uint8_t flat, recompflag;
  double sum_filtered, sum2_filtered, sum2, avgsum, stdsum, avg, std, slope, sigma_filter;


  /*  Loop through the temporary grid and filter the data.  */

  for (n = 0 ; n < height ; n++)
    {
      for (m = 0 ; m < width ; m++)
        {
          sumcount = 0;
          sum2 = 0.0;
          stdsum = 0.0;
          avgsum = 0.0;


          /*  Don't try to filter empty bins.  */

          if (grid[n][m].count)
            {
              /*  Get the information from the 8 cells surrounding this cell to compute the composite
                  standard deviation and average.  */

              for (i = n - 1 ; i <= n + 1 ; i++)
                {
                  for (j = m - 1 ; j <= m + 1 ; j++)
                    {
                      /*  Make sure each cell is in the area (edge effect).  */

                      if (i >= 0 && i < height && j >= 0 && j < width)
                        {
                          if (grid[i][j].count && !grid[i][j].cleared)
                            {
                              avgsum += grid[i][j].avg;
                              stdsum += grid[i][j].std;
                              sum2 += grid[i][j].avg * grid[i][j].avg;

                              sumcount++;
                            }
                        }
                    }
                }


              /*  Compute the eight slopes from the center cell to find out if it's flat enough to use the average
                  of the standard deviations or if we need to use the standard deviation of the averages.  We use a
                  reference slope of 1 degree to determine which we need to use.  */

              flat = NVTrue;
              for (i = n - 1 ; i <= n + 1 ; i++)
                {
                  for (j = m - 1 ; j <= m + 1 ; j++)
                    {
                      /*  Make sure each cell is in the area (edge effect).  */

                      if (i >= 0 && i < height && j >= 0 && j < width)
                        {
                          /*  Don't use the center cell.  */

                          if (i != n || j != m)
                            {
                              if (grid[i][j].count && ! grid[i][j].cleared)
                                {
                                  slope = (fabs (grid[n][m].avg - grid[i][j].avg)) / dx;
                      
                                  if (slope > 1.0)
                                    {
                                      flat = NVFalse;
                                      break;
                                    }
                                }
                            }
                        }
                    }
                }


              /*  If the slope is low (< 1 degree) we'll use an average of the cell standard deviations to beat the
                  depths against.  Otherwise, we'll compute the standard deviation from the cell averages.  Since we're
                  using the average of the computed standard deviations of all of the nine cells or the standard
                  deviations of the averages of all nine cells we multiply the resulting standard deviation (?) by two
                  to get a reasonable result, otherwise the standard deviation surface is too smooth and we end up
                  cutting out too much good data.  I must admit I arrived at these numbers by playing with the filter
                  using em3000 shallow water data and em121a deep water data but they appear to work properly.  This
                  way three sigma seems to cut out what you would expect three sigma to cut out.  If you leave it as is
                  it cuts out about 30%.  This is called empirically determining a value (From Nero's famous statement
                  "I'm the emperor and I can do what I damn well please, now hand me my fiddle.").   JCD  */

              avg = avgsum / (double) sumcount; 
              if (flat || sumcount < 2)
                {
                  std = (stdsum / (double) sumcount) * 2.0;
                }
              else
                {
                  std = (sqrt ((sum2 - ((double) sumcount * (avg * avg))) / ((double) sumcount - 1.0))) * 2.0;
                }

              sigma_filter = std_env * std;


              recompflag = NVFalse;
              for (i = 0 ; i < grid[n][m].count ; i++)
                {
                  grid[n][m].depths[i].filtered = NVFalse;


                  /*  Check for deep filter only.  */

                  if (deep)
                    {
                      if (adep[grid[n][m].depths[i].index] - avg >= sigma_filter) 
                        {
                          grid[n][m].depths[i].filtered = NVTrue;
                          recompflag = NVTrue;
                        }
                    }
                  else
                    {
                      if (fabs (adep[grid[n][m].depths[i].index] - avg) >= sigma_filter)
                        {
                          grid[n][m].depths[i].filtered = NVTrue;
                          recompflag = NVTrue;
                        }
                    }
                }


              if (recompflag)
                {
                  sum_filtered = 0.0;
                  sum2_filtered = 0.0;
                  filtered_count = 0;


                  for (i = 0 ; i < grid[n][m].count ; i++)
                    {
                      if (!(grid[n][m].depths[i].filtered))
                        {
                          sum_filtered += adep[grid[n][m].depths[i].index];
                          sum2_filtered += (adep[grid[n][m].depths[i].index] * adep[grid[n][m].depths[i].index]);

                          filtered_count++;
                        }
                    }


                  if (!filtered_count)
                    {
                      grid[n][m].cleared = NVTrue;
                    }
                  else
                    {
                      grid[n][m].avg = sum_filtered / (double) filtered_count; 
                      if (filtered_count > 1)
                        {
                          grid[n][m].std = sqrt ((sum2_filtered - ((double) filtered_count * 
                                                                   (pow ((double) grid[n][m].avg, 2.0)))) / 
                                                 ((double) filtered_count - 1.0));
                        }
                      else
                        {
                          grid[n][m].std = 0.0;
                        }
                    }
                }
            }
        }
    }
}
