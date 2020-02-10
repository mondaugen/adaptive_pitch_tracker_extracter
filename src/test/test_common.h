#ifndef COMMON_H
#define COMMON_H 

#include <stdlib.h>
#include <stdio.h>

#define PRINT_ANY(x)\
    printf(_Generic( (x), int : "%d", float : "%f"),x)

#define PRINT_ALL(a_,n) \
    ({ __auto_type a__ = (a_);\
       __auto_type n__ = (n);\
       while (n__--) { PRINT_ANY(*a__++); printf(" "); };\
       printf("\n"); })

#define CHK_EQ(a_,b_,n) \
    ({ __auto_type a__ = (a_);\
       __auto_type b__ = (b_);\
       __auto_type n__ = (n);\
       int eq__ = 1 ;\
       while (n__--) { eq__ &= (*a__++ == *b__++); };\
       eq__; })

struct err_msg {
    int err_num;
    char *err_name;
    struct err_msg *next;
};

static struct err_msg no_err_msg = {
    .err_num = 0,
    .err_name = "No error",
    .next = NULL,
};
static struct err_msg *err_msg_table = &no_err_msg;

static int gen_err_msg(char *msg)
{
    struct err_msg *old_err_msg_table = err_msg_table;
    err_msg_table = calloc(1,sizeof(struct err_msg));
    *err_msg_table = (struct err_msg) {
        .err_num = old_err_msg_table->err_num + 1,
        .err_name = msg,
        .next = old_err_msg_table
    };
    return err_msg_table->err_num;
}

static void print_err_msg(int err_num)
{
    struct err_msg *msg = err_msg_table;
    while (msg != NULL) {
        if (msg->err_num == err_num) {
            puts(msg->err_name);
            return;
        }
        msg = msg->next;
    }
    puts("Error number not found.");
}

static inline long
get_file_length(
    FILE *f)
{
    fseek(f,0,SEEK_END);
    long ret = ftell(f);
    rewind(f);
    return ret;
}

static inline long
get_file_length_path(const char *path)
{
    FILE *f = fopen(path,"r");
    if (!f) { return 0; }
    long length = get_file_length(f);
    fclose(f);
    return length;
}

static inline const char *
getenv_default(
    const char *env_name,
    const char *default_value)
{
    const char *ret = getenv(env_name);
    if (!ret) { return default_value; }
    return ret;
}

/*
if longest_length not 0, allocates that length, otherwise allocates length
that fits whole file
*/
static inline void *
file_to_array(const char *path, long longest_length)
{
    long file_length;
    void *file_contents = NULL;
    FILE *f = fopen(path,"r");
    if (!f) { goto fail; }
    if (longest_length != 0) {
        file_length = longest_length;
    } else {
        file_length = get_file_length(f);
    }
    file_contents = calloc(file_length,1);
    if (!file_contents) { goto fail; }
    fread(file_contents,1,file_length,f);
    fclose(f);
    return file_contents;
fail:
    if (f) { fclose(f); }
    if (file_contents) { free(file_contents); }
    return NULL;
}

static inline int
array_to_file(const char *path, const void *data, unsigned int length)
{
    FILE *f = fopen(path,"w");
    if (!f) { goto fail; }
    if ((length != fwrite(data,1,length,f))) { goto fail; };
    return 0;
fail:
    return -1;
}
    

#endif /* COMMON_H */
