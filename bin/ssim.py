#!/usr/bin/env python3 

import numpy as np
import cv2 as cv
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import peak_signal_noise_ratio as psnr
import sys


def calculate_ssim_psnr(image1, image2):
    """Calculates SSIM and PSNR between two images."""

    # Ensure images are in the same format (e.g., float64)
    image1 = image1.astype(np.uint8)
    image2 = image2.astype(np.uint8)

    # Calculate SSIM
    ssim_value = ssim(image1, image2, data_range=image1.max() - image1.min(), channel_axis=2)

    # Calculate PSNR
    psnr_value = psnr(image1, image2, data_range=image1.max() - image1.min())

    return ssim_value, psnr_value

# Example usage
if __name__ == "__main__":


    # Read images
    img1 = cv.imread(sys.argv[1])
    img2 = cv.imread(sys.argv[2])

    # Calculate SSIM and PSNR
    ssim_val, psnr_val = calculate_ssim_psnr(img1, img2)

    print(ssim_val)
    print(psnr_val)
