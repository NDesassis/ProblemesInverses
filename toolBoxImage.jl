# Filtre 
function prepare_filter(M,N,f0h,f0w,filtername)

    """
    
    M: image height
    N: image width
    f0h: cut-off frequency of the filter in the width axis
    f0w: cut-off frequency of the filter in the height axis
    filtername: "Gaussian" for periodized Gaussian, or "Cauchy" for Cauchy distribution in spectrum
    """
    

freqh=ifftshift((1:M).-(floor(M/2)+1)) # vecteur des fréquences (frequence 0 en 1)
freqw=ifftshift((1:N).-(floor(N/2)+1)) # vecteur des fréquences (frequence 0 en 1)

if filtername == "Gaussian"
    hath=[exp(-((fh/f0h)^2+(fw/f0w)^2)/2) for fh in freqh, fw in freqw] # "gaussienne" -> très régulier        
elseif filtername == "Cauchy"
    hath=[1 ./(1 .+(fh/f0h)^2+(fw/f0w)^2) for fh in freqh, fw in freqw] # Filtre en 1/k^2 -> continu mais non dérivable
    else
        error("Type de filter inconnu")
end
h=fftshift(real.(ifft(hath))) # Réponse impulsionnelle 
return h,hath,freqh,freqw
end
