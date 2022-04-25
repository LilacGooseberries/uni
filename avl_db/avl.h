#ifndef TREE_H
#define TREE_H

#include "student.h"

#define MAX_LVL 4

class tree_node:public student
{   
    private:
    
    tree_node *left;
    tree_node *right;
    
    tree_node *next;
    
    int bal;
    
    public:
        
    tree_node(const char *n = 0, int v = 0): student(n, v) { left = 0; right = 0; next = 0; bal = 0; }
    ~tree_node() { left = 0; right = 0; next = 0; bal = 0; }
        
    void print()
        { printf("%s %d (%d)\n", get_name(), get_value(), bal); }
        
    friend class tree;
};

class tree
{
    private:
    
    tree_node *root;
    
    public:
        
    tree() { root = 0; }
    ~tree();
    void destr(tree_node *r);
    
    tree_node *get_root()
        { return root; }
    
    int insert(tree_node *r, tree_node *add, tree_node *r_par = 0);
    
    int remove(tree_node *r, student *rem, tree_node *r_par = 0);
    
    int search(student *s);
    
    int balance(tree_node *r, tree_node **new_r);
    
    int L1(tree_node *a);
    int L2(tree_node *a);
    int R1(tree_node *a);
    int R2(tree_node *a);
    
    int insert_file(FILE *fp);
    void delete_file(FILE *fp);
    void search_file(FILE *fp);
    
    void print(tree_node *root, int lvl = 0);
};

#endif
