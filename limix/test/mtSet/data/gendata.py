import scipy as SP
import pdb
SP.random.seed(0)

if __name__=='__main__':

    N = 500
    P = 4
    S = 1000
    R = 10

    hr = 0.05
    hg = 0.40

    Xr = 1.*(SP.rand(N,R)<0.2)
    Xg = 1.*(SP.rand(N,S)<0.2)
    Xg-= Xg.mean(0)
    Xg/= Xg.std(0)
    XX = SP.dot(Xg,Xg.T)
    XX/= XX.diagonal().mean()
    XX+= 1e-4*SP.eye(XX.shape[0])

    Yr = SP.dot(Xr,SP.randn(R,P))
    Yr*= SP.sqrt(hr/Yr.var(0))
    Yg = SP.dot(Xg,SP.randn(S,P))
    Yg*= SP.sqrt(hg/Yg.var(0))
    Yn = SP.randn(N,P)
    Yn*= SP.sqrt((1-hr-hg)/Yn.var(0))

    Y = Yr+Yg+Yn
    Y-= Y.mean(0)
    Y/= Y.std(0)

    SP.savetxt('Xr.txt',Xr,fmt='%d')
    SP.savetxt('XX.txt',XX,fmt='%.6f')
    SP.savetxt('Y.txt',Y,fmt='%.6f')

