#include "avl.h"

int main()
{
    char buf[LEN], name[LEN], c;
    FILE *fp;
    tree a = tree();

    while(fgets(buf, LEN, stdin))
    {
        if(sscanf(buf, "%c %s", &c, name) != 2)
        {
            printf("incorrect command\n");
            break;
        }
        if(!(fp = fopen(name, "r")))
        { 
            printf("cannot open file\n");
            break; 
        }
        switch(c)
        {
            case 'i':
                
                a.insert_file(fp);
                break;
                
            case 'd':
                
                a.delete_file(fp);
                break;
                
            case 's':
                
                a.search_file(fp);
                break;
                
            default:
                
                printf("wrong command\n");
        }
        
        fclose(fp);
        /*printf("\n");
        a.print(a.get_root());
        printf("\n");*/
    }
    
    return 0;
}
