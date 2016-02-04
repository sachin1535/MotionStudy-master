function Hi = invH(H)
%INVH       Creates an inverse homogeneous transform for H

Hi = [H(1:3,1:3)',-H(1:3,1:3)'*H(1:3,4); 0,0,0,1];

end