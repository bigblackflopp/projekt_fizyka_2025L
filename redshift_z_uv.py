import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from PIL import Image

C = 299792.458  # speed of light in km/s
lista = []
lista_pr = []
def redshift(distance_mpc, H0=70):
    v = H0 * distance_mpc  # km/s
    beta = v / C
    if beta >= 1:
        return float("inf")
    z = np.sqrt((1 + beta) / (1 - beta)) - 1
    lista.append(distance_mpc)
    lista_pr.append(v)
    return z

def wavelength_to_rgb(wavelength):
    gamma = 0.8
    R = G = B = 0
    if 380 <= wavelength < 440:
        R = -(wavelength - 440) / (440 - 380)
        B = 1.0
    elif 440 <= wavelength < 490:
        G = (wavelength - 440) / (490 - 440)
        B = 1.0
    elif 490 <= wavelength < 510:
        G = 1.0
        B = -(wavelength - 510) / (510 - 490)
    elif 510 <= wavelength < 580:
        R = (wavelength - 510) / (580 - 510)
        G = 1.0
    elif 580 <= wavelength < 645:
        R = 1.0
        G = -(wavelength - 645) / (645 - 580)
    elif 645 <= wavelength <= 750:
        R = 1.0

    if 380 <= wavelength < 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (40)
    elif 420 <= wavelength < 645:
        factor = 1.0
    elif 645 <= wavelength <= 750:
        factor = 0.3 + 0.7 * (750 - wavelength) / (105)
    else:
        factor = 0.0

    R = pow(R * factor, gamma)
    G = pow(G * factor, gamma)
    B = pow(B * factor, gamma)

    return (int(R * 255), int(G * 255), int(B * 255))

def redshift_pixel(r, g, b, z):
    # Step 1: Simulate UV contribution for bright blue regions
    # Assume bright blue = UV-rich
    global uv_emission
    blue_strength = b - (r + g) / 2
    is_uv_source = blue_strength > 50 and b > 150 & r!=0

    # Assign synthetic UV wavelength
    if is_uv_source:
        uv_emission = 300  # nm (near UV)
    else:
        # Estimate visible emitted wavelength (basic)
        if r >= g and r >= b and r!=0 and g !=0 and b !=0:
            uv_emission = 620  # red
        elif g >= r and g >= b and r!=0 and g !=0 and b !=0:
            uv_emission = 530  # green
        elif b >=r and b >=g and r!=0 and g !=0 and b !=0:
            uv_emission = 460  # blue
        else:
            uv_emission = 0

    # Step 2: Redshift the wavelength
    shifted_wave = uv_emission * (1 + z)

    # Step 3: Convert to RGB or dim if outside visible range
    if shifted_wave > 750:
        # Far-infrared: nearly invisible
        return (int(r * 0.2), int(g * 0.1), int(b * 0.1))
    elif shifted_wave < 380:
        # Still in UV: invisible
        return (0, 0, 0)
    else:
        return wavelength_to_rgb(shifted_wave)


def adjust_brightness(img, z):
    arr = np.array(img).astype(np.float32) / 255.0  # normalize to [0,1]
    if z > 0.5:
        dim_factor = 1.0 - (z - 0.2)  # gradual dimming factor after z > 0.7
        arr *= dim_factor  # dim the entire image
        arr = np.clip(arr, 0, 1)  # Ensure values remain in [0, 1]
    return Image.fromarray((arr * 255).astype(np.uint8))

def redshift_image_spectral(img, z):
    img = img.convert("RGB")
    arr = np.array(img)
    out = np.zeros_like(arr)
    for y in range(arr.shape[0]):
        for x in range(arr.shape[1]):
            out[y, x] = redshift_pixel(*arr[y, x], z)
    img_out = Image.fromarray(out)
    return adjust_brightness(img_out, z)

def animate_redshift(img_path):
    original = Image.open(img_path)
    fig, ax = plt.subplots(figsize=(10,10))
    ax.axis('off')
    im = ax.imshow(original)

    #Andromedas distance to Earth is around 0.9 Mpc
    distances = np.linspace(1, 5000, 500)
    redshifts = [redshift(d) for d in distances]

    for d, z in zip(distances, redshifts):
        print(f"Distance: {d:.1f} Mpc → z = {z:.4f}")

    frames = []

    def update(frame):

        z = redshifts[frame]
        if (z<1.005) :
            global i
            i += 1
            shifted = redshift_image_spectral(original, z)
            im.set_data(shifted)
            ax.set_title(f"z = {z:.2f}, distance = {lista[i]}, galaxy 'speed' = {lista_pr[i]} km/s", fontsize=14)

    ani = FuncAnimation(fig, update, frames=len(redshifts), interval=10)
    plt.show()
    return distances, redshifts

def plot_redshift(distances, redshifts):
    plt.figure(figsize=(13, 10))
    plt.plot(distances, redshifts, 'r-')
    plt.xlabel("Distance (Mpc)")
    plt.ylabel("Redshift (z)")
    plt.ylim(0,2)
    plt.title("Redshift vs Distance (Relativistic)")
    plt.grid(True)
    plt.plot(distances,np.ones(500),"b-")
    plt.show()

if __name__ == "__main__":
    i = 0
    image_path = sys.argv[1]  # replace with your image
    distances, redshifts = animate_redshift(image_path)
    plot_redshift(distances, redshifts)
