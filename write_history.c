
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#ifdef NVWIN3X
    #include <winsock.h>
#endif

#include "gsf.h"
#include "nvutility.h"


int32_t write_history (int32_t argc, char **argv, char *comment, char *gsfFile __attribute__ ((unused)), int32_t handle)
{
    int             i;
    int             ret;
    int             len;
    extern int      gsfError;
    char            str[1024];
    time_t          t;
    gsfRecords      rec;
    gsfDataID       gsfID;


    memset (&rec, 0, sizeof(rec));


    /* Load the contents of the gsf history record */

    time (&t);
    rec.history.history_time.tv_sec = t;
    rec.history.history_time.tv_nsec = 0;

    str[0] = '\0';
    len = 0;
    for (i = 0 ; i < argc ; i++)
    {
        len += strlen (argv[i]) + 1;
        if ((unsigned int) len >= sizeof(str))
        {
            return(-1);
        }
        strcat (str, argv[i]);
        strcat (str, " ");
    }

    rec.history.command_line = str;

    gethostname (rec.history.host_name, sizeof (rec.history.host_name));

    rec.history.comment = comment;


    /* Seek to the end of the file and write a new history record */

    ret = gsfSeek (handle, GSF_END_OF_FILE);
    if (ret) fprintf(stderr, "gsfSeek error: %d\n", gsfError);

    memset (&gsfID, 0, sizeof(gsfID));
    gsfID.recordID = GSF_RECORD_HISTORY;
    ret = gsfWrite (handle, &gsfID, &rec);

    if (ret < 0)
    {
        return (gsfError);
    }


    return(0);
}
