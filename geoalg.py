
class GeoAlgBasis:
    def __init__(self,geoalg,name,ind):
        self.geoalg = geoalg
        self.name = name
        self.ind = ind

    def __str__(self):
        return self.name
    def __repr__(self):
        return str(self)

class GeoAlg:
    def __init__(self,dim,names=None,signs=None):
        if names == None:
            names = ['e'+str(i+1) for i in range(dim)]
        if signs == None:
            signs = [1 for i in range(dim)]
        self.signs = signs
        self.basis = [GeoAlgBasis(self,names[i],i) for i in range(dim)]
        self.blade = [GeoAlgBlade(self,1,[b]) for b in self.basis]
        self.multi = [GeoAlgMultiVector(self,[b]) for b in self.blade]
        self.scalar = GeoAlgMultiVector(self,[GeoAlgBlade(self,1,[])])
        self.I = GeoAlgMultiVector(self,[GeoAlgBlade(self,1,self.basis)])

    def getBasis(self):
        return self.multi

class GeoAlgBlade:
    def __init__(self,geoalg,coef,basis):
        self.geoalg = geoalg
        self.coef = coef
        self.basis = basis

    def copy(self):
        return GeoAlgBlade(self.geoalg,self.coef,self.basis.copy())

    def __eq__(self,other):
        s = self
        o = other
        if type(o)!=GeoAlgBlade:
            o = GeoAlgBlade(self.geoalg,o,[])
        if s.coef!=o.coef or len(s.basis)!=len(o.basis):
            return False
        for i in range(len(s.basis)):
            if s.basis[i]!=o.basis[i]:
                return False
        return True

    def __req__(self,other):
        s = self
        o = other
        if type(o)!=GeoAlgBlade:
            o = GeoAlgBlade(self.geoalg,o,[])
        if s.coef!=o.coef or len(s.basis)!=len(o.basis):
            return False
        for i in range(len(s.basis)):
            if s.basis[i]!=o.basis[i]:
                return False
        return True

    def __str__(self):
        s = str(self.coef)
        for b in self.basis:
            s+=str(b)
        return s
        
    def __repr__(self):
        return str(self)

class GeoAlgMultiVector:
    def __init__(self,geoalg,blades):
        self.geoalg = geoalg
        self.blades = blades

    def copy(self):
        blades = []
        for b in self.blades:
            blades.append(b.copy())
        return GeoAlgMultiVector(self.geoalg,blades)

    def simplify(self):
        blades = self.blades
        newblades = []
        for b in blades:
            basis = b.basis
            i = 0
            bas = len(basis)
            while i <= bas:
                j = 0
                while j < len(basis)-1:
                    if basis[j].ind>basis[j+1].ind:
                        basis[j],basis[j+1] = basis[j+1],basis[j]
                        b.coef = -b.coef
                    elif basis[j].ind == basis[j+1].ind:
                        b.coef *= b.geoalg.signs[basis[j].ind]
                        basis.pop(j)
                        basis.pop(j)
                    j+=1
                i+=1
        #print(self)
        i = 0
        while i < len(blades):
            #if blades[i].coef==0:
            #    blades.pop(i)
            #    continue
            j = i+1
            while j < len(blades):
                if blades[i].basis == blades[j].basis:
                    blades[i].coef += blades[j].coef
                    blades.pop(j)
                    continue
                j+=1
            if blades[i].coef==0:
                blades.pop(i)
                continue
            newblades.append(blades[i])
            i+=1
        if len(newblades)==0:
            newblades.append(GeoAlgBlade(self.geoalg,0,[]))
        self.blades = newblades

    def __neg__(self):
        mv = self.copy()
        for b in mv.blades:
            b.coef = -b.coef
        return mv

    def __add__(self,other):
        s = self.copy()
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        o = o.copy()
        s.blades.extend(o.blades)
        s.simplify()
        return s

    def __radd__(self,other):
        s = self.copy()
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        o = o.copy()
        o.blades.extend(s.blades)
        o.simplify()
        return o

    def __sub__(self,other):
        s = self.copy()
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        o = -o
        s.blades.extend(o.blades)
        s.simplify()
        return s

    def __rsub__(self,other):
        s = -self
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        o = o.copy()
        o.blades.extend(s.blades)
        o.simplify()
        return o

    def __or__(self,other):
        return ((self*other)+(other*self))/2

    def __ror__(self,other):
        return self|other

    def __xor__(self,other):
        return ((self*other)-(other*self))/2

    def __rxor__(self,other):
        return self^other
    
    def __mul__(self,other):
        s = self
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        blades = []
        for sb in s.blades:
            for ob in o.blades:
                b = sb.basis.copy()
                b.extend(ob.basis.copy())
                blades.append(GeoAlgBlade(self.geoalg,sb.coef*ob.coef,b))
        mv = GeoAlgMultiVector(self.geoalg,blades)
        mv.simplify()
        return mv

    def __rmul__(self,other):
        return self*other

    def inverse(self):
        s = self.copy()
        s2 = (s*s).blades[0].coef
        for b in s.blades:
            b.coef /= s2
        return s
    
    def __truediv__(self,other):
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        return self*o.inverse()

    def __rtruediv__(self,other):
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        return self.inverse()*o

    def __pow__(self,other):
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(geoalg,[o])
        return exp(ln(self)*o)

    def __rpow__(self,other):
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        return exp(ln(o)*self)

    def __lt__(self,other):
        s = self
        o = other
        I = self.geoalg.I
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        return (s^(o*I.inverse()))*I

    def __rlt__(self,other):
        s = self
        o = other
        I = self.geoalg.I
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        return (o^(s*I.inverse()))*I

    def __gt__(self,other):
        s = self
        o = other
        I = self.geoalg.I
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(geoalg,[o])
        return I*((I.inverse()*s)^o)

    def __rgt__(self,other):
        s = self
        o = other
        I = self.geoalg.I
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(geoalg,[o])
        return I*((I.inverse()*o)^s)

    def __eq__(self,other):
        s = self
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        if len(s.blades)!=len(o.blades):
            return False
        for i in range(len(s.blades)):
            if s.blades[i]!=o.blades[i]:
                return False
        return True

    def __req__(self,other):
        s = self
        o = other
        if type(o)==int or type(o)==float or type(o)==complex:
            o = GeoAlgBlade(self.geoalg,o,[])
            o = GeoAlgMultiVector(self.geoalg,[o])
        if len(s.blades)!=len(o.blades):
            return False
        for i in range(len(s.blades)):
            if s.blades[i]!=o.blades[i]:
                return False
        return True

    

    def __str__(self):
        s = ''
        if len(self.blades)==0:
            return '0'
        for i,blade in enumerate(self.blades):
            if i!=0:
                s+=' + '
            s+=str(blade)
        return s

    def __repr__(self):
        return str(self)

def exp(n,t=10):
    d = 1
    x = 1
    ex = 1
    for i in range(1,t+1):
        x *= n
        d *= i
        ex+=x/d
    return ex

def ln(n,t=10):
    s = (n-1)/(n+1)
    d = -1
    x = 1/s
    l=0
    for i in range(t):
        d+=2
        x*=s*s
        l+=x/d
    return 2*l

def ln2(n,t=10):
    ans = 0
    i = t
    while i>0:
        ans = (i*i*n)/(i+1-i*n+ans) 
        i-=1
    ans = n/(1+ans)
    return ans

def sin(n,t=10):
    s = n
    d = 1
    i2 = 1
    x = s
    ans = x
    for i in range(t):
        i2+=2
        d*=i2*(i2-1)
        x*=-s*s
        ans+=x/d
    return ans

def cos(n,t=10):
    s = n
    d = 1
    i2 = 0
    x = 1
    ans = 1
    for i in range(t):
        i2+=2
        d*=i2*(i2-1)
        x*=-s*s
        ans+=x/d
    return ans

if __name__=='__main__':
    geoalg = GeoAlg(3)
    e1,e2,e3 = geoalg.getBasis()
    s = geoalg.scalar
    geoalg2 = GeoAlg(4,names=['t','x','y','z'],signs=[1,-1,-1,-1])
    t,x,y,z = geoalg2.getBasis()
    geoalg3 = GeoAlg(1,names = ['e'],signs = [0])
    e = geoalg3.getBasis()[0]
    print(e1,e2,e3)
    print(t,x,y,z)
