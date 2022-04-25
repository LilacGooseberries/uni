#ifndef STUDENT_H
#define STUDENT_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define LEN 1234

class student
{
    private:
    
    char *name;
    int value;
    
    protected:
    
    int init(const char *n, int v);
    
    public:
    
    student(const char *n = 0, int v = 0);
    ~student();
    
    student & operator=(const student &x);
    
    int get_value() const
        { return value; }
    
    const char *get_name() const
        { return (const char *) name; }
    
    int operator<(const student &b) const
    {
        int c = strcmp(name, b.name);
        if(c < 0)
            return 1;
        else if(c > 0)
            return 0;
        else
            return value < b.value;
    }
        
    int operator==(const student &b) const
        { return (value == b.value) && (!strcmp(name, b.name)); }
    
    int read(FILE *fp);

    void print()
        { printf("%s %d\n", name, value); }
};

#endif
