#include "avl.h"

void tree::destr(tree_node *r)
{
    if(r->left)
        destr(r->left);
    if(r->right)
        destr(r->right);
    
    for(tree_node *p = r; p; )
    {
        tree_node *t = p->next;
        delete p;
        p = t;
    }
}

tree::~tree()
{
    if(root)
    {
        destr(root);
        root = 0;
    }
}

void tree::print(tree_node *root, int lvl)
{
    int i;
    
    if(!root || (lvl > MAX_LVL))
        return;
    
    for(i = 0; i < lvl; i++)
        printf("        ");
    
    root->print();
    
    print(root->left, lvl+1);
    print(root->right, lvl+1);
}

int tree::insert_file(FILE *fp)
{
    if(!root)
    {
        root = new tree_node();
        if(root->read(fp) == -1)
        {
            delete root;
            root = 0;
            return -1;
        }
        printf("inserted ");
        root->student::print();
    }
    
    while(1)
    { 
        tree_node *n = new tree_node();
        if(n->read(fp) == -1)
        {
            delete n;
            break;
        }
        
        insert(root, n);
        printf("inserted ");
        n->student::print();
    }

    return 0;
}

void tree::delete_file(FILE *fp)
{
    student *s = new student();
    while(s->read(fp) != -1)
    {
        if(remove(root, s) != -1)
        {
            printf("deleted ");
            s->print();
        }
    }
    delete s;
}

void tree::search_file(FILE *fp)
{
    student *s = new student();
    while(s->read(fp) != -1)
    {
        if(search(s) == 0)
        {
            printf("found ");
            s->print();
        }
    }
    delete s;
}

int tree::search(student *s)
{
    for(tree_node *p = root; p; )
    {
        if(*p < *s)      p = p->right;
        else if(*s < *p) p = p->left;
        else
        {
            tree_node *i;
            for(i = p; i && !(*i == *s); i = i->next);
            if(!i) return -1;
            return 0;
        }
    }
    return -1;
}

int tree::insert(tree_node *r, tree_node *add, tree_node *r_par)
{
    if(*add < *r)
    {
        if(!r->left)
        {
            r->left = add;
            (r->bal)--;
            return abs(r->bal);
        }
        else
        {
            if(!insert(r->left, add, r))
                return 0;
            
            (r->bal)--;
            switch(r->bal)
            {
                case 0: return 0;
                    
                case -1: return 1;
                    
                case -2:
                {
                    tree_node *new_root;
                    
                    int d = balance(r, &new_root);
                    
                    if(r_par)
                    {
                        if(*r < *r_par) r_par->left  = new_root;
                        else            r_par->right = new_root;
                    }
                    else
                        root = new_root;
                    
                    return d;
                }
                default: printf("dfq\n"); return -1;
            }
        }
    }
    else if(*r < *add)
    {
        if(!r->right)
        {
            r->right = add;
            (r->bal)++;
            return r->bal;
        }
        else
        {
            if(!insert(r->right, add, r))
                return 0;
            
            (r->bal)++;
            switch(r->bal)
            {
                case 0: return 0;
                    
                case 1: return 1;
                
                case 2:
                {
                    tree_node *new_root;
                    
                    int d = balance(r, &new_root);
                    
                    if(r_par)
                    {
                        if(*r < *r_par) r_par->left  = new_root;
                        else            r_par->right = new_root;
                    }
                    else
                        root = new_root;
                    
                    return d;
                }
                default: printf("dfq\n"); return -1;
            }
        }
    }
    else
    {
        tree_node *t = r->next;
        r->next = add;
        add->next = t;
        
        return 0;
    }
}

int tree::remove(tree_node *r, student *rem, tree_node *r_par)
{
    if(!r)
        return -1;
    if(*rem < *r)
    {
        int ret = remove(r->left, rem, r);
        switch(ret)
        {
            case -1: return -1;
                    
            case 0: return 0;
            
            case 1:
            {
                (r->bal)++;
                switch(r->bal)
                {
                    case 0: return 1;
                        
                    case 1: return 0;
                        
                    case 2:
                    {
                        tree_node *new_root;
                        
                        int b = balance(r, &new_root);
                        
                        if(r_par)
                        {
                            if(r_par->left == r) r_par->left  = new_root;
                            else                 r_par->right = new_root;
                        }
                        else
                            root = new_root;
                        
                        return 1-b;
                    }
                    default: printf("dfq\n"); return -1;
                }
            }
            default: printf("dfq\n"); return -1;
        }
    }
    else if(*r < *rem)
    {
        int ret = remove(r->right, rem, r);
        switch(ret)
        {
            case -1: return -1;
                    
            case 0: return 0;
            
            case 1:
            {
                (r->bal)--;
                switch(r->bal)
                {
                    case 0: return 1;
                        
                    case -1: return 0;
                        
                    case -2:
                    {
                        tree_node *new_root;
                        
                        int b = balance(r, &new_root);
                        
                        if(r_par)
                        {
                            if(r_par->left == r) r_par->left  = new_root;
                            else                 r_par->right = new_root;
                        }
                        else
                            root = new_root;
                        
                        return 1-b;
                    }
                    default: printf("dfq\n"); return -1;
                }
            }
            default: printf("dfq\n"); return -1;
        }
    }
    else
    {
        if(r->next)
        {
            if(*r == *rem)
            {
                if(r_par)
                {
                    if(r_par->left == r) r_par->left  = r->next;
                    else                 r_par->right = r->next;
                }
                else
                    root = r->next;
                
                r->next->left = r->left;
                r->next->right = r->right;
                
                delete r;
                return 0;
            }
            
            for(tree_node *p = r; p; p = p->next)
            {
                if(*p->next == *rem)
                {
                    p->next = p->next->next;
                    delete p->next;
                    return 0;
                }
            }
            return -1;
        }

        if(!(*r == *rem)) return -1;
        
        if(!r->left && !r->right)
        {
            if(r_par)
            {
                if(r_par->left == r) r_par->left  = 0;
                else                 r_par->right = 0;
            }
            else
                root = 0;

            delete r;
            return 1;
        }
        if(!r->left)
        {
            if(r_par)
            {
                if(r_par->left == r) r_par->left  = r->right;
                else                 r_par->right = r->right;
            }
            else
                root = r->right;
            
            delete r;
            return 1;
        }
        
        tree_node *sub, *sub_par = r;
        
        for(sub = r->left; sub->right; sub = sub->right)
            sub_par = sub;

        if(r_par)
        {
            if(r_par->left == r) r_par->left  = sub;
            else                 r_par->right = sub;
        }
        else
            root = sub;
        
        tree_node *s_l = sub->left;
        int s_b = sub->bal;
        
        sub->left = r->left;
        sub->right = r->right;
        sub->bal = r->bal;
        
        if(sub_par == r) sub->left = r;
        else    sub_par->right = r;
        
        r->left = s_l;
        r->right = 0;
        r->bal = s_b;

        int d_l = remove(sub->left, rem, sub);
        sub->bal += d_l;
        
        if(!d_l) return 0;
        
        switch(sub->bal)
        {
            case 0: return 1;
                
            case 1: return 0;
            
            case 2:
            {
                tree_node *new_root;
                
                int b = balance(sub, &new_root);
                
                if(r_par)
                {
                    if(r_par->left == sub) r_par->left  = new_root;
                    else                   r_par->right = new_root;
                }
                else
                    root = new_root;
                
                return 1-b;
            }
            default: printf("dfq\n"); return -1;
        }
    }
}

int tree::balance(tree_node *r, tree_node **new_r)
{
    if(r->bal == -2)
    {
        if(r->left->bal < 1)
        {
            *new_r = r->left;
            return L1(r);
        }
        
        *new_r = r->left->right;
        return L2(r);
    }

    if(r->right->bal > -1)
    {
        *new_r = r->right;
        return R1(r);
    }
    
    *new_r = r->right->left;
    return R2(r);
}

int tree::L1(tree_node *a)
{
    tree_node *l = a->left, *lr = l->right;
    
    l->right = a;
    a->left = lr;
    
    (l->bal)++;
    if(l->bal == 0)
    { 
        a->bal = 0;
        return 0;
    }

    a->bal = -1;
    return 1;
}

int tree::L2(tree_node *a)
{
    tree_node *l = a->left, *lr = l->right;
    tree_node *lrl = lr->left, *lrr = lr->right;
    
    lr->left = l; 
    lr->right = a;
    l->right = lrl;
    a->left = lrr;
    
    switch(lr->bal)
    {
        case -1:
        {
            l->bal = 0;
            a->bal = 1;
            break;
        }
        case 0:
        {
            l->bal = 0;
            a->bal = 0;
            break;
        }
        case 1:
        {
            l->bal = -1;
            a->bal = 0;
            break;
        }
    }
    lr->bal = 0;
    return 0;
}

int tree::R1(tree_node *a)
{
    tree_node *r = a->right, *rl = r->left;
    
    r->left = a;
    a->right = rl;
    
    (r->bal)--;
    if(r->bal == 0)
    { 
        a->bal = 0;
        return 0;
    }

    a->bal = 1;
    return 1;
}

int tree::R2(tree_node *a)
{
    tree_node *r = a->right, *rl = r->left;
    tree_node *rlr = rl->right, *rll = rl->left;
    
    rl->right = r; 
    rl->left = a;
    r->left = rlr;
    a->right = rll;
    
    switch(rl->bal)
    {
        case 1:
        {
            r->bal = 0;
            a->bal = -1;
            break;
        }
        case 0:
        {
            r->bal = 0;
            a->bal = 0;
            break;
        }
        case -1:
        {
            r->bal = 1;
            a->bal = 0;
            break;
        }
    }
    rl->bal = 0;
    return 0;
}
