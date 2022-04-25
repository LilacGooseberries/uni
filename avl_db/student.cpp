#include "student.h"

int student::init(const char *n, int v)
{
    value = v;
    if(n)
    {
        if(name)
        {
            delete [] name;
            name = 0;
        }
        name = new char [strlen(n)+1];
        if(!name)
            return -1;
        strcpy(name, n);
        return 0;
    }
    name = 0;
    return 0;
}

student::student(const char *n, int v)
{
    init(n, v);
}

student::~student()
{
    if(name)
    {
        delete [] name;
        name = 0;
    }
    value = 0;
}

student& student::operator=(const student &x)
{
    this->~student();
    init(x.name, x.value);
    return *this;
}

int student::read(FILE *fp)
{ 
    char buf[LEN];
    
    if(fscanf(fp, "%s %d", buf, &value) != 2)
    {
        if(feof(fp))
            return -1;
        return 1;
    }
    
    return init(buf, value);
}
