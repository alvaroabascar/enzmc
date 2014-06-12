function v = michaelis(S, p)
    Vm = p(1);
    Km = p(2);
    v = (Vm.*S)./(Km + S);
endfunction

function nlr(S, yi, guess)
    [v_fit_nl, PNL, cvg, iter] = leasqr(S', yi', guess', "michaelis");
    PNL
    iter
    plot(S, yi, 'or', v_fit_nl, '-g');
end

load Vmax
load Km
hist(Vmax, 300);
print('-dpng', 'Vmax_hist.png')
hist(Km, 300);
print('-dpng', 'Km_hist.png')
