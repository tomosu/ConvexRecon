#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "Phantom.h"

#define MBF 256

void delete_ws(char *string,
               char *ret_string)
{
    int i, j;

    j = 0;
    for (i = 0;; i++) {
        if (*(string + i) == '\0') {
            ret_string[j] = '\0';
            break;
        }
        if (  *(string + i) != '\t'
           && *(string + i) != ' '
           && *(string + i) != '\n') {
            ret_string[j] = (*(string + i));
            j++;
        }
    }
}

int replace_str(char *string,
                char *before,
                char *after)
{
    char  *find;
    size_t be_len = strlen(before);
    size_t af_len = strlen(after);

    if (be_len == 0 || (find = strstr(string, before)) == NULL) return 0;
    memmove(find + af_len, find + be_len, strlen(string) - (find + be_len - string) + 1);
    memcpy(find, after, af_len);
    return 1;
}

int count_body(char *filepath)
{
    FILE *fp = NULL;
    int   i;
    char  read_buf[MBF], ret_buf[MBF];

    if ((fp = fopen(filepath, "r")) == NULL) {
        printf("material file is missing\n");
        exit(0);
    }

    i = 0;
    while (fgets(read_buf, MBF, fp) != NULL) {
        delete_ws(read_buf, ret_buf);
        if (strstr(ret_buf, "BODY") != NULL) {
            i++;
        }
    }
    return i;
}

Phantom_t *setupPhantom(char *filepath)
{
    FILE      *fp          = NULL;
    Phantom_t *phantom_ary = NULL;

    int   body_count = count_body(filepath);
    char  read_buf[MBF], ret_buf[MBF];
    char *split_buf;

    phantom_ary = (Phantom_t *)malloc(sizeof(Phantom_t) * body_count);
    if ((fp = fopen(filepath, "r")) == NULL) {
        printf("phantom setup file is missing 1\n");
        exit(0);
    }

    int j = 0;
    while (fgets(read_buf, MBF, fp) != NULL) {
        delete_ws(read_buf, ret_buf);

        if (strstr(ret_buf, "BODY") != NULL) {
            j++;
        }

        if (j < 1) {
            printf("phantom setup file is invalid\n");
            exit(0);
        } else  {
            if (strstr(ret_buf, "RAD=") != NULL) {
                replace_str(ret_buf, "RAD=", "");
                phantom_ary[j - 1].rad = atof(ret_buf);
            }

            if (strstr(ret_buf, "MU=") != NULL) {
                replace_str(ret_buf, "MU=", "");
                phantom_ary[j - 1].mu = atof(ret_buf);
            }

            if (strstr(ret_buf, "POSITION=") != NULL) {
                replace_str(ret_buf, "POSITION=", "");
                replace_str(ret_buf, "(", "");
                replace_str(ret_buf, ")", "");

                Real2_t v;
                double     tmp;
                split_buf = strtok(ret_buf, ",");
                tmp       = atof(split_buf);
                v.x       = tmp;
                split_buf = strtok(NULL, ",");
                tmp       = atof(split_buf);
                v.y       = tmp;

                phantom_ary[j - 1].center = v;
            }

            if (strstr(ret_buf, "IS_BASE=") != NULL) {
                replace_str(ret_buf, "IS_BASE=", "");
                phantom_ary[j - 1].isBase = atoi(ret_buf);
            }
        }
    }

    return phantom_ary;
}

#ifdef _TEST_SETUP_PHANTOM_
int main()
{
    char phantomPath[128] = "phantoms/phantom.ini";

    Phantom_t *phantoms;
    phantoms = setupPhantom(phantomPath);

    for (int i = 0; i < 4; i++) {
        printf("rad:%f\n", phantoms[i].rad);
        printf("rad:%f\n", phantoms[i].mu);
        printf("pos:%f %f\n", phantoms[i].center.x, phantoms[i].center.y);
        printf("base:%d\n", phantoms[i].isBase);
    }
    return 0;
}

#endif // ifdef _TEST_SETUP_PHANTOM_
