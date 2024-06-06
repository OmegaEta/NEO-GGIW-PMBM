function shape_=predictShapeRotate(pre_v,post_v,shape)
    if(norm(pre_v)==0||norm(post_v)==0)
        shape_=shape;
        return;
    end
    d=2;
    
    angle=acos(dot(pre_v,post_v)/(norm(pre_v)*norm(post_v)));

    if(det([pre_v post_v])<0)
        angle=-angle;
    end
    
    R=[cos(angle) -sin(angle);sin(angle) cos(angle)];
    T=[R zeros(2,2);zeros(2,2) R];
    
    shape_=shape;
    
    for i=1:length(shape_)
        shape_(i).m=T*shape_(i).m;
        shape_(i).P=T*shape_(i).P*T';
        
        S=shape_(i).V/(shape_(i).v-2*d-2);
        S = (S + S')/2;
        
        shape_(i).V=(shape_(i).v-2*d-2)*R*S*R';
    end
end