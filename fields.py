import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from joblib import Parallel, delayed
import scipy
import os
import pandas as pd


class CustomPDF(scipy.stats.rv_continuous):
    def __init__(self, h_t, h, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.h_t = h_t  # Store h_t inside self
        self.h = h      # Store h inside self

    def _pdf(self, d):
        mu1 = 4.7849
        sigma1 = 0.7022
        term1 = (1 / (sigma1 * np.sqrt(2 * np.pi))) * np.exp(-((np.log(d) - mu1) ** 2) / (2 * sigma1 ** 2))
        mu2 = 1.3386
        sigma2 = 0.9418
        term2 = (1 / (sigma2 * np.sqrt(2 * np.pi))) * np.exp(-((np.log(d) - mu2) ** 2) / (2 * sigma2 ** 2))
        pdf = term1 * (1 - Heaviside(self.h, self.h_t)) + term2 * Heaviside(self.h, self.h_t)
        return pdf / d
    
    def _rvs(self, size, random_state=None):
        rng = self._random_state if random_state is None else random_state

        # Ensure size is an integer
        size = size[0]

        batch_size = max(10000, size * 10)

        samples = np.empty([size])
        
        max_pdf = 100
        envelope = max_pdf*(1/(self.b - self.a))

        accepted = 0
        while accepted<size:

            y = np.random.uniform(self.a,self.b,size = batch_size)

            pdf_y = self._pdf(y)

            u = np.random.rand(batch_size)

            accept = u < pdf_y / envelope
            if np.any(accept):
                accepted_candidates = y[accept]
                num_accept = accepted_candidates.size
                remaining = size - accepted
                if num_accept > remaining:
                    accepted_candidates = accepted_candidates[:remaining]
                samples[accepted:accepted + accepted_candidates.size] = accepted_candidates
                accepted += accepted_candidates.size

        return samples
    
class superparticle:
    def __init__(self):
        # define constants
        self.mu_0 = 1.2566370621219e-6
        self.eps_0 = 8.854187812813e-12
        self.q_e = 1.60217663e-19
        
        # other particle properties
        self.r = np.NaN #grain radius
        self.q = np.NaN #grain charge
        self.R = np.NaN #rotation radius
        self.h = np.NaN #rotation height
        self.Vt = np.NaN #tangential velocity of rotation
        self.w = [0, 0, np.NaN] #angular velocity of rotation
        self.x = np.NaN #x position
        self.y = np.NaN #y position
    
    def fieldAtPoint(self, point):
        """calculate electric and magnetic field of particle at a specific point"""
        c = np.array([self.x, self.y, self.h])
        distance = -c + point

        B = (self.mu_0 * self.q_e * self.q) / (4 * np.pi * np.linalg.norm(distance)**3) * np.cross(np.cross(self.w, c), distance)
        E = (self.q_e * self.q) / (4 * np.pi * self.eps_0 * np.linalg.norm(distance)**3) * distance

        return B, E
    
    def particletonp(self):
        """output particle characteristics as array"""
        return np.array([self.r, self.R, self.q, self.h, self.Vt, self.x, self.y])
    
    def timestep(self, t, v):
        """update particle characteristics based on timestep t and dust devil velocity v"""
        rotation = self.w[2] * t 
        #new position = old position + rotation + distance travelled
        self.x = self.x + (self.R * np.cos(rotation)) + (v[0] * t)
        self.y = self.y + (self.R * np.sin(rotation)) + (v[1] * t)
    
def concat(a, b, c):
    """Concatenates three arrays column-wise into a single 2D array."""
    return np.concatenate([np.reshape(a,[-1,1]), np.reshape(b,[-1,1]), np.reshape(c,[-1,1])],1)

def rad2height(D):
    """Calculates the height based on the given dust devil radius."""
    return ((D/2) * (1/0.913)) **2

def dd_volume(D, H):
    """Computes the volume of a cylindrical dust devil."""
    outer_radius = D / 2
    #assuming all dust devil inside as Metzgar
    inner_radius = D / 4
    return H *  np.pi * (outer_radius**2 - inner_radius**2)

def grain_diameter(h,nparts):
    """sample from the grain diameter distribution based on height."""
    h_t = 2.0 

    custom_dist = CustomPDF(h_t, h, a=0.0001, b=1000, name="custom")

    samples = custom_dist.rvs(size = nparts)

    return samples

def dd_rad_profile(D):
    """Generates a random radius within the dust devil profile."""
    radius = (D/2) * ((1) - (0.5 * (np.random.rand())))
    return radius

def tan_vel_profile(R):
    """calculates tangential velocity based on motion radius"""
    vel = 0.0007982*R**2 + -0.1637*R + 10.68
    return vel

def charge_profile(radius):
    """parameterised charge based on radius"""
    #ed 1
    #q = 12.53*radius**2 - 3.73e4 *radius**-0.8

    #ed 0.5
    #q = 6.272*radius**2 - 2.027e4 *radius**-0.8

    #ed 0.45
    #q = 5.643*radius**2 - 1.936e4 *radius**-0.8

    #ed 2
    #q = 25.07*radius**2 - 5.464e4 *radius **-0.8

    #ed 2.5
    # q = 31.33*radius**2 - 6.76e4 *radius **-0.8

    #ed 3
    q = 37.6*radius**2 - 1.039e5 *radius **-0.8

    #ed 10
    #q = 125.3*radius**2 - 1.585e5 *radius **-0.8

    return q

def Heaviside(h,a):
    """smoothed heaviside switch function"""
    return 0.5 * np.tanh(h-a)+0.5
    
def initialise(D, N_Sim, N_AtHeight):
    """initialise each particle into a list of custom class superparticles"""
    dust_rho = 2.0e6  # 2.0g/cm3 to g/m3

    devil_centre = np.array([0, 0, 0])

    #calculate max height and thus dust loaded volume
    max_height = rad2height(D)
    volume = dd_volume(D, max_height)

    #calculate mass of dust and number of particle assuming average particle radius of 30 micron
    dust_loading = 0.296  # mg/m3 to g/m3
    dust_mass = dust_loading * volume
    av_mass = dust_rho * 4/3 * np.pi * (3e-6) ** 3

    #determine average number of particles in dust devil
    N = np.ceil(dust_mass / av_mass)

    clump_factor = N / N_Sim

    N_heights = int(N_Sim / N_AtHeight)

    heights = np.linspace(0,max_height,N_heights)

    particles = [superparticle() for _ in range(N_Sim)]
    with tqdm(total=N_heights * N_AtHeight, desc="Initializing Particles", leave=False) as init_bar:
        for i in range(N_heights):
            diameters = grain_diameter(heights[i],N_AtHeight)
            for j in range(N_AtHeight):
                idx = (i * N_AtHeight) + j
                particles[idx].h = heights[i]
                particles[idx].r = diameters[j]/2
                particles[idx].R = dd_rad_profile(D)
                particles[idx].Vt = tan_vel_profile(particles[idx].R)
                particles[idx].w[2] = particles[idx].Vt / particles[idx].R * -1
                theta_offset = 2*np.pi*np.random.rand()
                particles[idx].x = devil_centre[0] + particles[idx].R * np.cos(theta_offset)
                particles[idx].y = devil_centre[1] +   particles[idx].R * np.sin(theta_offset)
                particles[idx].q = clump_factor * charge_profile(particles[idx].r) * -1

                init_bar.update(1)

    return particles

def compute_field(q, pt):
    return q.fieldAtPoint(pt)

def parallelSim(FWHM, N_Sim, N_AtHeight):
    """Parallelized simulation of field computation for multiple particles at one observation point."""
    
    p = initialise(FWHM, N_Sim, N_AtHeight)
    pt = np.array([0, 0, 0])

    B_total = np.array([0, 0, 0])
    E_total = np.array([0, 0, 0])

    results = Parallel(n_jobs=4, backend="loky")(
        delayed(compute_field)(q, pt) for q in p
    )

    fields = np.array(results)

    # Compute total fields
    B_total = np.sum(fields[:,0],axis=0)
    E_total = np.sum(fields[:,1],axis=0)

    return B_total, E_total

def simfunc(FWHM, N_Sim, N_AtHeight):
    """one core sim at points, good for 1e4 1e5 particles, okay for 1e6, slow for 1e7"""
    p = initialise(FWHM, N_Sim, N_AtHeight)

    points = [np.array([0, 0, 0])]

    for pt in tqdm(points, desc="Points", leave=False):
        fields = []
        for q in tqdm(p, desc="Particles", leave=False):
            [Bpart, Epart] = q.fieldAtPoint(pt)
            fields.append([Bpart, Epart])
    
    fields = np.array(fields)

def farrell_event(FWHM, N_Sim, N_AtHeight, point):
    """one core sim at point"""
    p = initialise(FWHM, N_Sim, N_AtHeight)
    
    fields = []
    for q in tqdm(p, desc="Particles", leave=False):
        [Bpart, Epart] = q.fieldAtPoint(point)
        fields.append([Bpart, Epart])
        
    fields = np.array(fields)

    B_total = np.sum(fields[:,0],axis=0)
    E_total = np.sum(fields[:,1],axis=0)
    return B_total, E_total, p

def saveparticles(parts):
    data = [q.particletonp() for q in parts]
    index = pd.Index(range(len(data)), dtype='int64')
    output = pd.DataFrame(data, columns=['r', 'R', 'q', 'h', 'Vt', 'x', 'y'], index=index)
    output.to_csv('particles.csv')

def compute_field_at_points(p, pt):
    """"Vectorized computation for multiple points at once."""
    positions = np.array([[q.x, q.y, q.h] for q in p])
    charges = np.array([q.q for q in p])
    w_vectors = np.array([q.w for q in p])

    mu_0 = p[0].mu_0
    q_e = p[0].q_e
    eps_0 = p[0].eps_0

    distances = pt[:, np.newaxis, :] - positions[np.newaxis, :, :]
    norms = np.linalg.norm(distances, axis=2)

    cross_w_c = np.cross(np.cross(w_vectors[np.newaxis, :, :], positions[np.newaxis, :, :]), distances)
    B = (mu_0 * q_e * charges[np.newaxis, :, np.newaxis]) / (4 * np.pi * norms[:, :, np.newaxis]**3) * cross_w_c

    E = (q_e * charges[np.newaxis, :, np.newaxis]) / (4 * np.pi * eps_0 * norms[:, :, np.newaxis]**3) * distances

    B_total = np.sum(B, axis=1)  # (M, 3)
    E_total = np.sum(E, axis=1)  # (M, 3)

    return np.hstack((pt, B_total, E_total))




def main():
    #np.random.seed(42)

    FWHM = 7 #metres
    N_AtHeight = 100
    
    ###convergence/sensitivity simulations
    ###convergence
    # N_Sim = [int(1e3), int(1e4), int(1e5), int(1e6)]
    # n_runs = 1000
    # ###sensitivity
    # # N_Sim = [int(1e6)]
    # # n_runs = 100

    # for m_Sim in N_Sim:
    #     output = []
    #     its = 0
    #     with tqdm(total=n_runs, desc="Simulation Runs") as pbar:
    #         while its < n_runs:
    #             #parallel
    #             output.append([m_Sim, parallelSim(FWHM, m_Sim, N_AtHeight)])
    #             #single core
    #             #output.append([m_Sim, simfunc(FWHM, m_Sim, N_AtHeight)])
    #             its += 1
    #             pbar.update(1)
    #     df = pd.DataFrame([a[0], a[1][0][0], a[1][0][1], a[1][0][2], a[1][1][0], a[1][1][1], a[1][1][2]] for a in output)
    #     df.columns = ["N", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]    
    #     df.to_csv(f'convergence_sims_{m_Sim}.csv')
    #     df.to_csv(f'sensitivity_-dl.csv')


    ###farrell verification sim
    # n_runs = 100
    # N_Sim = int(1e6)
    # output_centre = []
    # output_50m = []
    # its = 0
    # with tqdm(total=n_runs, desc="Simulation Runs") as pbar:
    #     while its < n_runs:
    #         output_centre.append([N_Sim, farrell_event(FWHM, N_Sim, N_AtHeight, np.array([0,0,1]))])
    #         output_50m.append([N_Sim, farrell_event(FWHM, N_Sim, N_AtHeight, np.array([50,0,1]))])
    #         its += 1
    #         pbar.update(1)
    # df = pd.DataFrame([a[0], a[1][0][0], a[1][0][1], a[1][0][2], a[1][1][0], a[1][1][1], a[1][1][2]] for a in output_centre)
    # df.columns = ["N", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]
    # df.to_csv('centre.csv')
    # df = pd.DataFrame([a[0], a[1][0][0], a[1][0][1], a[1][0][2], a[1][1][0], a[1][1][1], a[1][1][2]] for a in output_50m)
    # df.columns = ["N", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]
    # df.to_csv('50m.csv')

    ###contour map 
    # N_Sim = int(1e6)
    # p = initialise(FWHM, N_Sim, N_AtHeight)
    # saveparticles(p)
    # xs = np.linspace(-50,50,1001)
    # #ys = np.array([0])
    # ys = np.linspace(-50,49.9,1000)
    # #zs = np.linspace(0,49,500)
    # zs = np.array([1])
    # X, Y, Z = np.meshgrid(xs,ys,zs)
    # points = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T
    # batch_size = len(ys)//10
    # num_batches = len(points) // batch_size + 1
    # for i in tqdm(range(num_batches), desc="Processing Batches"):
    #     output = compute_field_at_points(p,points[i * batch_size:(i + 1) * batch_size])
    #     batch_filename = f"contour_part_{i}.csv"
    #     df = pd.DataFrame(output)
    #     df.columns = ["x", "y", "z", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]
    #     df.to_csv(batch_filename)

    ###design sweep
    # FWHM = np.linspace(0.2,1.0,17)
    # n_runs = 1
    # N_Sim = int(1e6)
    # output_centre = []
    # for size in tqdm(FWHM, desc="Sizes"):
    #     its = 0
    #     with tqdm(total=n_runs, desc="Simulation Runs") as pbar:
    #         while its < n_runs:
    #             output_centre.append([size, farrell_event(size, N_Sim, N_AtHeight, np.array([0,0,0]))])
    #             its += 1
    #             pbar.update(1)
    # df = pd.DataFrame([a[0], a[1][0][0], a[1][0][1], a[1][0][2], a[1][1][0], a[1][1][1], a[1][1][2]] for a in output_centre)
    # df.columns = ["D", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]
    # df.to_csv('sizing_cases_0m.csv')


    FWHM = np.array([7])
    n_runs = 1
    N_Sim = int(1e6)
    output_centre = []
    for size in tqdm(FWHM, desc="Sizes"):
        its = 0
        with tqdm(total=n_runs, desc="Simulation Runs") as pbar:
            while its < n_runs:
                output_centre.append([size, farrell_event(size, N_Sim, N_AtHeight, np.array([50,0,1]))])
                its += 1
                pbar.update(1)
    df = pd.DataFrame([a[0], a[1][0][0], a[1][0][1], a[1][0][2], a[1][1][0], a[1][1][1], a[1][1][2]] for a in output_centre)
    df.columns = ["D", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]
    p  = output_centre[0][1][2]
    
    #df.to_csv('sizing_cases_0m.csv')

if __name__ == "__main__":
    main()
    
plt.show()