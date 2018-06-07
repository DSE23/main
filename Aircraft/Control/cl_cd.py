def nonlinear(aoa):
    contin = True
    continCD = True
    f = open("cl_alpha.txt")
    cl_alpha = f.readlines()
    f.close
    
    g = open("cd_alpha.txt")
    cd_alpha = g.readlines()
    g.close
    
    cl_alpha_lst = len(cl_alpha)*[0.]
    cl_lst = len(cl_alpha)*[0.]
    
    cd_alpha_lst = len(cd_alpha)*[0.]
    cd_lst = len(cd_alpha) * [0.]
    
    for i in range (0,len(cl_alpha)):
        cl_alpha[i] = cl_alpha[i].strip("\n")
        variables  = cl_alpha[i].split("\t")
        cl_alpha_lst[i] = variables[0]
        cl_lst[i]       = variables[1]
        
    for i in range (0,len(cd_alpha)):
        cd_alpha[i] = cd_alpha[i].strip("\n")
        variables  = cd_alpha[i].split("\t")
        cd_alpha_lst[i] = variables[0]
        cd_lst[i]       = variables[1]

    for i in range (0,len(cl_alpha)):
        if float(cl_alpha_lst[i]) > float(aoa) and contin:            
            dy_dx = (float(cl_lst[i])-float(cl_lst[i-1]))/(float(cl_alpha_lst[i])-float(cl_alpha_lst[i-1]))
        
            newcl = float(cl_lst[i-1]) + (float(aoa) - float(cl_alpha_lst[i-1]))* dy_dx
            contin = False
            
    for i in range (0,len(cd_alpha)):
        if float(cd_alpha_lst[i]) > float(aoa) and continCD:            
            dy_dx = (float(cd_lst[i])-float(cd_lst[i-1]))/(float(cd_alpha_lst[i])-float(cd_alpha_lst[i-1]))
        
            newcd =  float(cd_lst[i-1]) + (float(aoa) - float(cd_alpha_lst[i-1]))* dy_dx
            continCD = False
    return(newcl,newcd)
    

    