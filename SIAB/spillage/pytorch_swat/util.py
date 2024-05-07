import torch
import torch_complex

def ND_list(*sizes,element=None):
    size_1,*size_other = sizes
    l = [element] * size_1
    if size_other:
        for i in range(len(l)):
            l[i] = ND_list(*size_other,element=element)
    else:
        if element in ["dict()","list()"]:    
            for i in range(size_1):    
                l[i] = eval(element)
    return l

def ndlist(*dims, element = None):
    pass    

def ignore_line(file,N):
    for _ in range(N):    
        file.readline()
        
        
class Info:
    def Nm(self,il): return 2*il+1
    def __str__(self):
        return "\n".join([name+"\t"+str(value) for name,value in self.__dict__.items()])
    __repr__=__str__

def convert_tocuda(src):
    """convert a casscatte of data to cuda, for list, dict, torch.Tensor, torch_complex.ComplexTensor"""
    if isinstance(src, list):
        return [convert_tocuda(x) for x in src]
    if isinstance(src, dict):
        return {i: convert_tocuda(x) for i, x in src.items()}
    if isinstance(src, torch.Tensor):
        return src.cuda()
    if isinstance(src, torch_complex.ComplexTensor):
        return torch_complex.ComplexTensor(src.real.cuda(), src.imag.cuda())
    raise TypeError(f"Unexpected type of input: {type(src)}, only support list, dict, torch.Tensor, torch_complex.ComplexTensor and their nested structure.")

def update0(t):
    return t.masked_fill(mask=(t==0), value=1E-10)

def Nm(il):
    return 2*il+1