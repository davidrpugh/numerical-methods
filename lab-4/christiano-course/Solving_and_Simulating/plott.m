load C:\sw20\docs\ECB\paper2005\AER_revision\johannsen\cmrfiles\Second_order\cee\linear outp infl cons inv
outplin=outp;
infllin=infl;
conslin=cons;
invlin=inv;

load C:\sw20\docs\ECB\paper2005\AER_revision\johannsen\cmrfiles\Second_order\cee\nopruning outp infl cons inv
outpnoprune=outp;
inflnoprune=infl;
consnoprune=cons;
invnoprune=inv;

load C:\sw20\docs\ECB\paper2005\AER_revision\johannsen\cmrfiles\Second_order\cee\pruning outp infl cons inv
outpprune=outp;
inflprune=infl;
consprune=cons;
invprune=inv;

tt=1:length(outp);

subplot(221)
plot(tt,outplin,tt,outpnoprune,'*-',tt,outpprune,'o-')
legend('linear','no pruning','pruning')
title('output')
axis tight

subplot(222)
plot(tt,infllin,tt,inflnoprune,'*-',tt,inflprune,'o-')
title('inflation')
axis tight

subplot(223)
plot(tt,conslin,tt,consnoprune,'*-',tt,consprune,'o-')
title('consumption')
axis tight

subplot(224)
plot(tt,invlin,tt,invnoprune,'*-',tt,invprune,'o-')
title('investment')
axis tight
