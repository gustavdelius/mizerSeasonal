# Gonadic mass dynamics

## Deriving the PDE for gonadic mass

Standard mizer describes the number density \\N(w,t)\\ by the PDE \\
\frac{\partial N}{\partial t} = -\frac{\partial}{\partial
w}(g\\N)-\mu\\N, \tag{1}\\ where \\g(w,t)\\ is the growth rate of an
individual of size \\w\\ and \\\mu(w,t)\\ is its death rate. This
equation simply states that the rate of change of the number of
individuals in an infinitesimal size range around \\w\\ is given by the
difference between the rate at which individuals grow into the size
range and the rate at which they grow out (the transport term), minus
the rate at which they die (the loss term).

Instead of the number density \\N\\ we could have used the biomass
density \\B = w\\N\\ to describe the population. Its evolution equation
can easily be derived from [Equation 1](#eq-PDEN) \\\begin{split}
\frac{\partial B}{\partial t} &= w\frac{\partial N}{\partial t} =
-w\frac{\partial}{\partial w}\left(g\\\frac{B}{w}\right)-\mu\\B,\\
&=-\frac{\partial}{\partial w}\left(g\\B\right)+g\\\frac{B}{w}-\mu\\B.
\end{split}\\ This is again intuitive: in addition to the transport term
and the loss term we now also have a source term \\g B/w = gN\\
representing coming from the growth taking place within the
infinitesimal size class.

If we now denote by \\Q(w,t)\\ the total gonadic biomass density and by
\\G(w,t)\\ the gonadic growth rate of an individual of size \\w\\ (when
it is not reproducing), then they satisfy the analogous equation \\
\frac{\partial Q}{\partial t} =-\frac{\partial}{\partial
w}\left(g\\Q\right)+G\\N-\mu\\Q-r\\Q, \\ where we have an additional
loss term representing the release of gonadic mass during reproduction.
The release rate \\r(w,t)\\ will encode the seasonality of this
reproduction.

In mizer, the total energy \\E_r\\ available for growth an reproduction
is split as \\g=E_r(1-\psi)\\ and \\G=E_r\psi\\. Thus
\\G=g\\\psi/(1-\psi)\\.

Instead of considering the total gonadic biomass density \\Q(w,t)\\, we
could consider the gonadic mass \\q(w,t)\\ of an individual fish of size
\\w\\, given by \\ q= \frac{Q}{N}. \\ Let’s derive its evolution PDE:
\\\begin{split} \frac{\partial q}{\partial t} &=
\frac{1}{N}\frac{\partial Q}{\partial t}- \frac{Q}{N^2}\frac{\partial
N}{\partial t}\\ &=\frac{1}{N}\left(-\frac{\partial}{\partial
w}\left(g\\Q\right)+G\\N-\mu\\Q-r\\Q\right)
-\frac{Q}{N^2}\left(-\frac{\partial}{\partial w}(g\\N)-\mu\\N\right).
\end{split} \tag{2}\\ To simplify this we observe that the two mortality
terms cancel and that \\\begin{split}
\frac{Q}{N^2}\frac{\partial}{\partial w}(g\\N)
&=\frac{q}{N}\frac{\partial}{\partial w}\left(g\\\frac{Q}{q}\right)\\
&=\frac{1}{N}\frac{\partial}{\partial w}\left(g\\Q\right)
-g\\\frac{\partial q}{\partial w}. \end{split}\\ Substituting this into
[Equation 2](#eq-PDEwG) and cancelling terms gives \\ \frac{\partial
q}{\partial t}=G-g\\\frac{\partial q}{\partial w}-r\\ q. \tag{3}\\ At
this point we feel a bit stupid when we realise that \\
\frac{dq(w,t)}{dt}=\frac{\partial q(w,t)}{\partial t}+ \frac{\partial
q(w,t)}{\partial w}\frac{dw}{dt} \\ and that \\dw/dt = g\\, so that
[Equation 3](#eq-wg) reduces to \\ \frac{dq}{dt}=G-r\\q, \tag{4}\\ which
is what we could have written down immediately, given our definition of
\\G\\. So this equation would have been a good starting point for
deriving the PDE [Equation 3](#eq-wg).

All the above equations are species specific. So in particular there is
one gonadic mass function \\q_i(w,t)\\ for each species \\i\\ and we
will need to specify one rate function \\r_i(w,t)\\ for each species
\\i\\.

We need to specify the mass-specific release rate \\r(w,t)\\. This then
determines the total reproduction rate \\\begin{equation}
R_p(t)=\frac{\epsilon}{2w_0}\int N(w,t)q(w,t)r(w,t)dw. \end{equation}\\
It will be difficult to know how to choose the mass-specific rate
\\r(w,t)\\ to produce the observed egg production rate \\R_p(t)\\. We’ll
probably have to do this by trial and error initially. The rate
\\r(w,t)\\ should become large at least once a year so that the gonadic
mass goes back to close to zero. One possibility would be to make
\\r(w,t)\\ independent of \\w\\.

## Numerical implementation

To solve the PDE [Equation 3](#eq-wg) we use the same upwind-difference
scheme used by mizer to solve the PDE [Equation 1](#eq-PDEN). Let us
therefore recall the details of that scheme.

We discretise [Equation 1](#eq-PDEN) as follows: \\
\frac{N_w^{t+1}-N_w^t}{\Delta
t}=-\frac{g_w^tN_w^{t+1}-g\_{w-1}^tN\_{w-1}^{t+1}}{\Delta
w}-\mu_w^t\\N_w^{t+1}. \tag{5}\\ Here we use \\w\\ and \\t\\ to denote
the integer indices of the discretised weight and discretised time
respectively. Multiplying by \\\Delta t\\ and collecting terms gives \\
N_w^{t+1}\left(1+g_w^t\frac{\Delta t}{\Delta w}+\mu_w^t\\\Delta
t\right)= N\_{w-1}^{t+1}\left(g\_{w-1}^t\frac{\Delta t}{\Delta w}\right)
+N_w^t. \tag{6}\\ Solving for \\N_w^{t+1}\\ gives \\
N_w^{t+1}=\frac{S_w^t-A_w^tN\_{w-1}^{t+1}}{B_w^t}, \tag{7}\\ where \\
\begin{split} A_w^t&=-g\_{w-1}^t\frac{\Delta t}{\Delta w},\\
B_w^t&=1+g_w^t\frac{\Delta t}{\Delta w}+\mu^t_w\\\Delta t,\\
S_w^t&=N_w^t. \end{split} \tag{8}\\ Given an initial value and a
boundary condition encoding the influx of new-born individuals to the
smallest size class, The recurrence relation [Equation 7](#eq-nds)
allows mizer to calculate \\N_w^{t+1}\\ for all size classes at all
later times.

Now we use the same scheme for the gonadic dynamics. So we discretise
[Equation 3](#eq-wg) as follows: \\ \frac{q_w^{t+1}-q_w^t}{\Delta
t}=-g_w^t\frac{q_w^{t+1}-q\_{w-1}^{t+1}}{\Delta
w}-r_w^t\\q_w^{t+1}+G_w^t. \\ Note how close this equation for \\q\\ is
to equation [Equation 5](#eq-nd) for \\N\\. Multiplying by \\\Delta t\\
and collecting terms gives \\ q_w^{t+1}\left(1+g_w^t\frac{\Delta
t}{\Delta w}+r_w^t\\\Delta t\right)=
q\_{w-1}^{t+1}\left(g_w^t\frac{\Delta t}{\Delta w}\right)
+q_w^t+G_w^t\\\Delta t. \\ Solving for \\q_w^{t+1}\\ gives \\
q_w^{t+1}=\frac{s_w^t-a_w^tq\_{w-1}^{t+1}}{b_w^t}, \tag{9}\\ where \\
\begin{split} a_w^t&=-g_w^t\frac{\Delta t}{\Delta w},\\
b_w^t&=1+g_w^t\frac{\Delta t}{\Delta w}+r_w^t\\\Delta t,\\
s_w^t&=q_w^t+G_w^t\\\Delta t. \end{split} \tag{10}\\ Together with an
initial value and the boundary condition that \\q_1^t=0\\ at all times
because the smallest individuals do not yet have any gonadic mass, this
recurrence relation allows us to calculate \\q_w^{t+1}\\ for all size
classes at all later times.

Note how small the difference is between the discretised equation
[Equation 9](#eq-qds) for \\q\\ and the equation [Equation 7](#eq-nds)
for \\N\\. The only difference lies in replacing \\A_w^t, B_w^t\\ and
\\S_w^t\\ given in [Equation 8](#eq-ABS) by \\a_w^t, b_w^t\\ and
\\s_w^t\\ given in [Equation 10](#eq-abs).

## Open problems

1.  Once we are modelling cohorts in mizer, the diffusion terms in the
    size-spectrum dynamics will become important, because they will lead
    to the broadening of the cohorts.

2.  We do not yet know what the required step sizes are in size and in
    time to faithfully reproduce the dynamics. Too large size steps will
    lead to too much numerical diffusion. Too large time steps may lead
    to numerical instability.
