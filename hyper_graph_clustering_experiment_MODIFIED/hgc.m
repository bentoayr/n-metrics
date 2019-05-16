function [p,ng_eff,cluster] = hgc(n,m,msize,mlist,ng,NR)
%
%	Alexei Vazquez, April 8, 2009
%	
%	Reference material:
%		Vazquez A Finding hypergraph communities: 
%		Bayesian approach and variational solution.
%		J. Stat. Mech. (2009) P07006
%	
%
%	hgc: hypergraph clustering function
%
%		hgc finds the subdivision of a hypergraph into clusters
%
%	Hypergraph Input
%
%		n	number of vertices
%		m	number of edges (hyper-edges)
%		msize	one dimensional array
%				msize(j) number of elements in edge j
%				msize(1), ..., mszie(m)
%
%		mlist	two-dimensional array
%				mlist(i) list containing the vertices belonging to edge i 
%				mlist(1,1), ..., mlist(1,msize(1))
%
%				mlist(j,1), ..., mlist(j,msize(j))
%
%				mlist(m,1), ..., mlist(m,msize(m))
%
%		Note (graph clustering): there are n edges, and edge i represents the list of nearest neighbors of a vertex i
%
%
%	Parameters Input
%
%		ng	maximum number of clusters
%		NR	number of sampling initial conditions
%
%	Output
%
%		p	two-dimensional array
%			probability p(i,k) that vertex i belongs to group k
%
%		ng_eff	inferred number of clusters
%
%		cluster	one-dimensional array assigning samples to clusters (after removal of empty clusters)
%
%	Warning! This code has not been tested extensively.
%

	% Set algorithm parameters

	rerror=1e-6;		% Relative error stopping threshold
	tilde_alpha=1e-6;	% prior exponent (non-informative limit, or Jayne's prior, tilde_alpha->0)
	tilde_beta=1e-6;	% prior exponent (non-informative limit, or Jayne's prior, tilde_beta->0)
	tilde_gamma=1e-6;	% prior exponent (non-informative limit, or Jayne's prior, tilde_gamma->0)
	tilde_alpha_plus_beta=tilde_alpha+tilde_beta;

	% Create nlist (list of edges in which each vertex participates, the reverse of mlist)

	nsize=zeros(n);
	for j=1:m,
		for s=1:msize(j),
			i=mlist(j,s);
			nsize(i)=nsize(i)+1;
			nlist(i,nsize(i))=j;
		end
	end

	% Save state of the random number generator

	state=rand('twister');

	% Recursive clustering algorithm, sampling NR initial conditions. Selects the solution with lowest F.

	Fmin=0;
	for r=1:NR,
		[p1,F1]=hgc1(n,nsize,nlist,ng,m,msize,mlist,n*ng*r);
		if r==1 || F1<Fmin,
			Fmin=F1;
			p=p1;
		end
	end

	% Assign best group

	group_index=zeros(ng);
	cluster=zeros(n);

	ng_eff=0;
	for i=1:n,
		best=1;
		for k=2:ng,
			if p(i,k)>p(i,best),
				best=k;
			end
		end
		if group_index(best)==0,
			ng_eff=ng_eff+1;
			group_index(best)=ng_eff;
		end
		cluster(i)=group_index(best);
	end

	% Restore the state of the random number generator

	rand('twister',state);

	% Function declarations

	function psi = digamma(x1,x2),
	% Series expasion of psi(x1)-psi(x2)
		psi=0;
		for i=1:100;
			psi=psi+(x1-x2)/((x1+i-1)*(x2+i-1));
		end
	end	% psi

	function [p,F] = hgc1(n,nsize,nlist,ng,m,msize,mlist,seed),
	%
	% Recursive clustering algorithm, sampling one initial condition
	%

		% Initialize p

		rand('twister',seed);

		for i=1:n,
			sum=0;
			for k=1:ng,
				p(i,k)=rand;
				sum=sum+p(i,k);
			end
			for k=1:ng, 
				p(i,k)=p(i,k)/sum;
			end
		end

		F0=1;
		F=2*F0;

		sgamma=ng*tilde_gamma+n;

		% iterations

		while 2*abs(F-F0)>(abs(F)+abs(F0))*rerror,

			% Update gamma and <log pi>

			psi0=psi(sgamma);
			for k=1:ng,
				gamma_(k)=tilde_gamma;
				for i=1:n, 
					gamma_(k)=gamma_(k)+p(i,k);
				end
				logpi(k)=digamma(gamma_(k),sgamma);
			end

			% Update alpha, beta, <log theta> and <log (1-theta)>

			for k=1:ng,
				psi0=psi(gamma_(k));
				for j=1:m,
					alpha(k,j)=tilde_alpha;
					for s=1:msize(j),
						i=mlist(j,s);
						alpha(k,j)=alpha(k,j)+p(i,k);
					end
					alpha_plus_beta=tilde_alpha_plus_beta+(gamma_(k)-tilde_gamma);
					logtheta(k,j)=digamma(alpha(k,j),alpha_plus_beta);
					log1_theta(k,j)=digamma(alpha_plus_beta-alpha(k,j),alpha_plus_beta);
				end
			end

			% Update p

			for i=1:n,
				for k=1:ng,
					h(k)=logpi(k);
					for j=1:m,
						h(k)=h(k)+log1_theta(k,j);
					end
					for s=1:nsize(i),
						j=nlist(i,s);
						h(k)=h(k)+logtheta(k,j)-log1_theta(k,j);
					end
					if k==1 || h(k)>hmax, 
						hmax=h(k); 
					end
				end
				sum=0;
				for k=1:ng,
					sum=sum+exp(h(k)-hmax);
				end
				for k=1:ng,
					p(i,k)=exp(h(k)-hmax)/sum;
				end
			end

			% Update optimization objetive

			F0=F;
			F=0;

			for i=1:n,
				for k=1:ng,
					if p(i,k)>0,
						F=F+p(i,k)*log(p(i,k));
					end
				end
			end

			for k=1:ng,
				for j=1:m,
					alpha_plus_beta=tilde_alpha_plus_beta+(gamma_(k)-tilde_gamma);
					F=F+gammaln(alpha_plus_beta)-gammaln(alpha(k,j))-gammaln(alpha_plus_beta-alpha(k,j));
				end
			end

			F=F+gammaln(sgamma);
			for k=1:ng,
				F=F-gammaln(gamma_(k));
			end

		end

	end	% hgc1

end	% hgc
