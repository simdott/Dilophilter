/*
 * Dilophilter LV2 Plugin - Dual filter plugin
 * Copyright (C) 2026 Simon Delaruotte
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */

#include <lv2/core/lv2.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdlib.h>
#include <stdio.h>

#define DILOPHILTER_URI "urn:simdott:dilophilter"

typedef enum {
    INPUT_LEFT = 0,
    INPUT_RIGHT = 1,
    OUTPUT_LEFT = 2,
    OUTPUT_RIGHT = 3,
    HP_SLOPE = 4,
    HP_CUTOFF = 5,
    HP_RESONANCE = 6,
    LP_SLOPE = 7,
    LP_CUTOFF = 8,
    LP_RESONANCE = 9,
} PortIndex;

typedef struct {
    // Ports
    const float* input_left;
    const float* input_right;
    const float* hp_slope;
    const float* hp_cutoff;
    const float* hp_resonance;
    const float* lp_slope;
    const float* lp_cutoff;
    const float* lp_resonance;
    float* output_left;
    float* output_right;
    
    // High-pass filters
    float hp12_a1, hp12_a2, hp12_b0, hp12_b1, hp12_b2;
    float hp12_x1[2], hp12_x2[2], hp12_y1[2], hp12_y2[2];
    
    float hp24_a1[2], hp24_a2[2], hp24_b0[2], hp24_b1[2], hp24_b2[2];
    float hp24_x1[2][2], hp24_x2[2][2], hp24_y1[2][2], hp24_y2[2][2];
    
    float hp48_a1[4], hp48_a2[4], hp48_b0[4], hp48_b1[4], hp48_b2[4];
    float hp48_x1[2][4], hp48_x2[2][4], hp48_y1[2][4], hp48_y2[2][4];
    
    float hp24_bw_a1[2], hp24_bw_a2[2], hp24_bw_b0[2], hp24_bw_b1[2], hp24_bw_b2[2];
    float hp24_bw_x1[2][2], hp24_bw_x2[2][2], hp24_bw_y1[2][2], hp24_bw_y2[2][2];
    
    float hp48_bw_a1[4], hp48_bw_a2[4], hp48_bw_b0[4], hp48_bw_b1[4], hp48_bw_b2[4];
    float hp48_bw_x1[2][4], hp48_bw_x2[2][4], hp48_bw_y1[2][4], hp48_bw_y2[2][4];
    
    // Low-pass filters
    float lp12_a1, lp12_a2, lp12_b0, lp12_b1, lp12_b2;
    float lp12_x1[2], lp12_x2[2], lp12_y1[2], lp12_y2[2];
    
    float lp24_a1[2], lp24_a2[2], lp24_b0[2], lp24_b1[2], lp24_b2[2];
    float lp24_x1[2][2], lp24_x2[2][2], lp24_y1[2][2], lp24_y2[2][2];
    
    float lp48_a1[4], lp48_a2[4], lp48_b0[4], lp48_b1[4], lp48_b2[4];
    float lp48_x1[2][4], lp48_x2[2][4], lp48_y1[2][4], lp48_y2[2][4];
    
    float lp24_bw_a1[2], lp24_bw_a2[2], lp24_bw_b0[2], lp24_bw_b1[2], lp24_bw_b2[2];
    float lp24_bw_x1[2][2], lp24_bw_x2[2][2], lp24_bw_y1[2][2], lp24_bw_y2[2][2];
    
    float lp48_bw_a1[4], lp48_bw_a2[4], lp48_bw_b0[4], lp48_bw_b1[4], lp48_bw_b2[4];
    float lp48_bw_x1[2][4], lp48_bw_x2[2][4], lp48_bw_y1[2][4], lp48_bw_y2[2][4];
    
    float prev_hp_cutoff, prev_hp_resonance;
    float prev_lp_cutoff, prev_lp_resonance;
    float sample_rate;
    float prev_hp_slope, prev_lp_slope;

} Dilophilter;

// High-pass Q conversion functions
static float hp_resonance_to_q_12dB(float resonance) {
    return resonance;
}

static float hp_resonance_to_q_24dB(float resonance) {
    const float butterworth_q = 0.5412f;
    const float target_q = 0.707f;
    const float scale = target_q / butterworth_q;
    return butterworth_q * scale * (resonance / 0.707f);
}

static float hp_resonance_to_q_48dB(float resonance, int section) {
    const float butterworth_q[4] = {0.5098f, 0.6013f, 0.8999f, 2.5628f};
    const float target_q = 0.707f;
    const float scale = target_q / butterworth_q[section];
    return butterworth_q[section] * scale * (resonance / 0.707f);
}

// Low-pass Q conversion functions
static float lp_resonance_to_q_12dB(float resonance) {
    return resonance;
}

static float lp_resonance_to_q_24dB(float resonance) {
    const float butterworth_q = 0.5412f;
    const float target_q = 0.707f;
    const float scale = target_q / butterworth_q;
    return butterworth_q * scale * (resonance / 0.707f);
}

static float lp_resonance_to_q_48dB(float resonance, int section) {
    const float butterworth_q[4] = {0.5098f, 0.6013f, 0.8999f, 2.5628f};
    const float target_q = 0.707f;
    const float scale = target_q / butterworth_q[section];
    return butterworth_q[section] * scale * (resonance / 0.707f);
}

// High-pass coefficient calculations
static void calculate_hp12_coefficients(Dilophilter* self, float cutoff, float resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;
    
    float Q = hp_resonance_to_q_12dB(resonance);
    if (Q < 0.1f) Q = 0.1f;
    if (Q > 10.0f) Q = 10.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    const float alpha = sin_omega / (2.0f * Q);
    const float a0 = 1.0f + alpha;
    
    const float b0 = (1.0f + cos_omega) / 2.0f;
    const float b1 = -(1.0f + cos_omega);
    const float b2 = (1.0f + cos_omega) / 2.0f;
    const float a1 = -2.0f * cos_omega;
    const float a2 = 1.0f - alpha;

    self->hp12_b0 = b0 / a0;
    self->hp12_b1 = b1 / a0;
    self->hp12_b2 = b2 / a0;
    self->hp12_a1 = a1 / a0;
    self->hp12_a2 = a2 / a0;
}

static void calculate_hp24_coefficients(Dilophilter* self, float cutoff, float resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);

    for (int i = 0; i < 2; ++i) {
        float Q = hp_resonance_to_q_24dB(resonance);
        if (Q < 0.1f) Q = 0.1f;
        if (Q > 10.0f) Q = 10.0f;
        
        const float alpha = sin_omega / (2.0f * Q);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f + cos_omega) / 2.0f;
        const float b1 = -(1.0f + cos_omega);
        const float b2 = (1.0f + cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->hp24_b0[i] = b0 / a0;
        self->hp24_b1[i] = b1 / a0;
        self->hp24_b2[i] = b2 / a0;
        self->hp24_a1[i] = a1 / a0;
        self->hp24_a2[i] = a2 / a0;
    }
}

static void calculate_hp48_coefficients(Dilophilter* self, float cutoff, float resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);

    for (int i = 0; i < 4; ++i) {
        float Q = hp_resonance_to_q_48dB(resonance, i);
        if (Q < 0.1f) Q = 0.1f;
        if (Q > 10.0f) Q = 10.0f;
        
        const float alpha = sin_omega / (2.0f * Q);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f + cos_omega) / 2.0f;
        const float b1 = -(1.0f + cos_omega);
        const float b2 = (1.0f + cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->hp48_b0[i] = b0 / a0;
        self->hp48_b1[i] = b1 / a0;
        self->hp48_b2[i] = b2 / a0;
        self->hp48_a1[i] = a1 / a0;
        self->hp48_a2[i] = a2 / a0;
    }
}

static void calculate_hp24_bw_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    const float Q_values[2] = {0.5412f, 1.3066f};

    for (int i = 0; i < 2; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f + cos_omega) / 2.0f;
        const float b1 = -(1.0f + cos_omega);
        const float b2 = (1.0f + cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->hp24_bw_b0[i] = b0 / a0;
        self->hp24_bw_b1[i] = b1 / a0;
        self->hp24_bw_b2[i] = b2 / a0;
        self->hp24_bw_a1[i] = a1 / a0;
        self->hp24_bw_a2[i] = a2 / a0;
    }
}

static void calculate_hp48_bw_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    const float Q_values[4] = {0.5098f, 0.6013f, 0.8999f, 2.5628f};

    for (int i = 0; i < 4; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f + cos_omega) / 2.0f;
        const float b1 = -(1.0f + cos_omega);
        const float b2 = (1.0f + cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->hp48_bw_b0[i] = b0 / a0;
        self->hp48_bw_b1[i] = b1 / a0;
        self->hp48_bw_b2[i] = b2 / a0;
        self->hp48_bw_a1[i] = a1 / a0;
        self->hp48_bw_a2[i] = a2 / a0;
    }
}

// Low-pass coefficient calculations
static void calculate_lp12_coefficients(Dilophilter* self, float cutoff, float resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;
    
    float Q = lp_resonance_to_q_12dB(resonance);
    if (Q < 0.1f) Q = 0.1f;
    if (Q > 10.0f) Q = 10.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    const float alpha = sin_omega / (2.0f * Q);
    const float a0 = 1.0f + alpha;
    
    const float b0 = (1.0f - cos_omega) / 2.0f;
    const float b1 = 1.0f - cos_omega;
    const float b2 = (1.0f - cos_omega) / 2.0f;
    const float a1 = -2.0f * cos_omega;
    const float a2 = 1.0f - alpha;

    self->lp12_b0 = b0 / a0;
    self->lp12_b1 = b1 / a0;
    self->lp12_b2 = b2 / a0;
    self->lp12_a1 = a1 / a0;
    self->lp12_a2 = a2 / a0;
}

static void calculate_lp24_coefficients(Dilophilter* self, float cutoff, float resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);

    for (int i = 0; i < 2; ++i) {
        float Q = lp_resonance_to_q_24dB(resonance);
        if (Q < 0.1f) Q = 0.1f;
        if (Q > 10.0f) Q = 10.0f;
        
        const float alpha = sin_omega / (2.0f * Q);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f - cos_omega) / 2.0f;
        const float b1 = 1.0f - cos_omega;
        const float b2 = (1.0f - cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->lp24_b0[i] = b0 / a0;
        self->lp24_b1[i] = b1 / a0;
        self->lp24_b2[i] = b2 / a0;
        self->lp24_a1[i] = a1 / a0;
        self->lp24_a2[i] = a2 / a0;
    }
}

static void calculate_lp48_coefficients(Dilophilter* self, float cutoff, float resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);

    for (int i = 0; i < 4; ++i) {
        float Q = lp_resonance_to_q_48dB(resonance, i);
        if (Q < 0.1f) Q = 0.1f;
        if (Q > 10.0f) Q = 10.0f;
        
        const float alpha = sin_omega / (2.0f * Q);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f - cos_omega) / 2.0f;
        const float b1 = 1.0f - cos_omega;
        const float b2 = (1.0f - cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->lp48_b0[i] = b0 / a0;
        self->lp48_b1[i] = b1 / a0;
        self->lp48_b2[i] = b2 / a0;
        self->lp48_a1[i] = a1 / a0;
        self->lp48_a2[i] = a2 / a0;
    }
}

static void calculate_lp24_bw_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    const float Q_values[2] = {0.5412f, 1.3066f};

    for (int i = 0; i < 2; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f - cos_omega) / 2.0f;
        const float b1 = 1.0f - cos_omega;
        const float b2 = (1.0f - cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->lp24_bw_b0[i] = b0 / a0;
        self->lp24_bw_b1[i] = b1 / a0;
        self->lp24_bw_b2[i] = b2 / a0;
        self->lp24_bw_a1[i] = a1 / a0;
        self->lp24_bw_a2[i] = a2 / a0;
    }
}

static void calculate_lp48_bw_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    const float Q_values[4] = {0.5098f, 0.6013f, 0.8999f, 2.5628f};

    for (int i = 0; i < 4; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        
        const float b0 = (1.0f - cos_omega) / 2.0f;
        const float b1 = 1.0f - cos_omega;
        const float b2 = (1.0f - cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->lp48_bw_b0[i] = b0 / a0;
        self->lp48_bw_b1[i] = b1 / a0;
        self->lp48_bw_b2[i] = b2 / a0;
        self->lp48_bw_a1[i] = a1 / a0;
        self->lp48_bw_a2[i] = a2 / a0;
    }
}

static LV2_Handle instantiate(
    const LV2_Descriptor* descriptor,
    double sample_rate,
    const char* bundle_path,
    const LV2_Feature* const* features) {
    
    Dilophilter* self = (Dilophilter*)malloc(sizeof(Dilophilter));
    if (!self) return NULL;
    
    // Initialize all high-pass filter states to zero
    for (int ch = 0; ch < 2; ++ch) {
        self->hp12_x1[ch] = self->hp12_x2[ch] = self->hp12_y1[ch] = self->hp12_y2[ch] = 0.0f;
    }
    
    for (int i = 0; i < 2; ++i) {
        for (int ch = 0; ch < 2; ++ch) {
            self->hp24_x1[ch][i] = self->hp24_x2[ch][i] = self->hp24_y1[ch][i] = self->hp24_y2[ch][i] = 0.0f;
            self->hp24_bw_x1[ch][i] = self->hp24_bw_x2[ch][i] = self->hp24_bw_y1[ch][i] = self->hp24_bw_y2[ch][i] = 0.0f;
        }
    }
    
    for (int i = 0; i < 4; ++i) {
        for (int ch = 0; ch < 2; ++ch) {
            self->hp48_x1[ch][i] = self->hp48_x2[ch][i] = self->hp48_y1[ch][i] = self->hp48_y2[ch][i] = 0.0f;
            self->hp48_bw_x1[ch][i] = self->hp48_bw_x2[ch][i] = self->hp48_bw_y1[ch][i] = self->hp48_bw_y2[ch][i] = 0.0f;
        }
    }
    
    // Initialize all low-pass filter states to zero
    for (int ch = 0; ch < 2; ++ch) {
        self->lp12_x1[ch] = self->lp12_x2[ch] = self->lp12_y1[ch] = self->lp12_y2[ch] = 0.0f;
    }
    
    for (int i = 0; i < 2; ++i) {
        for (int ch = 0; ch < 2; ++ch) {
            self->lp24_x1[ch][i] = self->lp24_x2[ch][i] = self->lp24_y1[ch][i] = self->lp24_y2[ch][i] = 0.0f;
            self->lp24_bw_x1[ch][i] = self->lp24_bw_x2[ch][i] = self->lp24_bw_y1[ch][i] = self->lp24_bw_y2[ch][i] = 0.0f;
        }
    }
    
    for (int i = 0; i < 4; ++i) {
        for (int ch = 0; ch < 2; ++ch) {
            self->lp48_x1[ch][i] = self->lp48_x2[ch][i] = self->lp48_y1[ch][i] = self->lp48_y2[ch][i] = 0.0f;
            self->lp48_bw_x1[ch][i] = self->lp48_bw_x2[ch][i] = self->lp48_bw_y1[ch][i] = self->lp48_bw_y2[ch][i] = 0.0f;
        }
    }
    
    self->prev_hp_slope = 5.0f;  // Default to "Off"
    self->prev_lp_slope = 5.0f;  // Default to "Off"

    self->sample_rate = (float)sample_rate;
    self->prev_hp_cutoff = 0.0f;
    self->prev_hp_resonance = 0.707f;
    self->prev_lp_cutoff = 0.0f;
    self->prev_lp_resonance = 0.707f;
    
    // Initialize all filters with default values
    calculate_hp12_coefficients(self, 30.0f, 0.707f);
    calculate_hp24_coefficients(self, 30.0f, 0.707f);
    calculate_hp48_coefficients(self, 30.0f, 0.707f);
    calculate_hp24_bw_coefficients(self, 30.0f);
    calculate_hp48_bw_coefficients(self, 30.0f);
    
    calculate_lp12_coefficients(self, 10000.0f, 0.707f);
    calculate_lp24_coefficients(self, 10000.0f, 0.707f);
    calculate_lp48_coefficients(self, 10000.0f, 0.707f);
    calculate_lp24_bw_coefficients(self, 10000.0f);
    calculate_lp48_bw_coefficients(self, 10000.0f);
    
    return (LV2_Handle)self;
}

static void connect_port(LV2_Handle instance, uint32_t port, void* data) {
    Dilophilter* self = (Dilophilter*)instance;
    
    switch ((PortIndex)port) {
        case INPUT_LEFT: 
            self->input_left = (const float*)data; 
            break;
        case INPUT_RIGHT: 
            self->input_right = (const float*)data; 
            break;
        case OUTPUT_LEFT: 
            self->output_left = (float*)data; 
            break;
        case OUTPUT_RIGHT: 
            self->output_right = (float*)data; 
            break;
        case HP_SLOPE: 
            self->hp_slope = (const float*)data; 
            break;
        case HP_CUTOFF: 
            self->hp_cutoff = (const float*)data; 
            break;
        case HP_RESONANCE: 
            self->hp_resonance = (const float*)data; 
            break;
        case LP_SLOPE: 
            self->lp_slope = (const float*)data; 
            break;
        case LP_CUTOFF: 
            self->lp_cutoff = (const float*)data; 
            break;
        case LP_RESONANCE: 
            self->lp_resonance = (const float*)data; 
            break;
    }
}

static void run(LV2_Handle instance, uint32_t n_samples) {
    Dilophilter* self = (Dilophilter*)instance;
    
    int hp_slope_selector = (int)*(self->hp_slope);
    int lp_slope_selector = (int)*(self->lp_slope);
    float current_hp_cutoff = *(self->hp_cutoff);
    float current_hp_resonance = *(self->hp_resonance);
    float current_lp_cutoff = *(self->lp_cutoff);
    float current_lp_resonance = *(self->lp_resonance);
        
    // Check if any parameters changed
    int parameters_changed = (current_hp_cutoff != self->prev_hp_cutoff) ||
                            (current_lp_cutoff != self->prev_lp_cutoff) ||
                            (current_hp_resonance != self->prev_hp_resonance) ||
                            (current_lp_resonance != self->prev_lp_resonance) ||
                            (hp_slope_selector != (int)self->prev_hp_slope) ||
                            (lp_slope_selector != (int)self->prev_lp_slope);
    
    if (hp_slope_selector == 5 && lp_slope_selector == 5 && !parameters_changed) {
    for (uint32_t i = 0; i < n_samples; i++) {
        self->output_left[i] = self->input_left[i];
        self->output_right[i] = self->input_right[i];
    }
    return;
}
    // Update coefficients only if parameters changed AND filter is not off
    if (parameters_changed) {
        // Update HPF coefficients if needed (skip if HPF is off)
        if (hp_slope_selector != 5) {
            switch (hp_slope_selector) {
                case 0: calculate_hp12_coefficients(self, current_hp_cutoff, current_hp_resonance); break;
                case 1: calculate_hp24_coefficients(self, current_hp_cutoff, current_hp_resonance); break;
                case 2: calculate_hp48_coefficients(self, current_hp_cutoff, current_hp_resonance); break;
                case 3: calculate_hp24_bw_coefficients(self, current_hp_cutoff); break;
                case 4: calculate_hp48_bw_coefficients(self, current_hp_cutoff); break;
            }
        }
        
        // Update LPF coefficients if needed (skip if LPF is off)
        if (lp_slope_selector != 5) {
            switch (lp_slope_selector) {
                case 0: calculate_lp12_coefficients(self, current_lp_cutoff, current_lp_resonance); break;
                case 1: calculate_lp24_coefficients(self, current_lp_cutoff, current_lp_resonance); break;
                case 2: calculate_lp48_coefficients(self, current_lp_cutoff, current_lp_resonance); break;
                case 3: calculate_lp24_bw_coefficients(self, current_lp_cutoff); break;
                case 4: calculate_lp48_bw_coefficients(self, current_lp_cutoff); break;
            }
        }
        
        // Update all previous values
        self->prev_hp_cutoff = current_hp_cutoff;
        self->prev_hp_resonance = current_hp_resonance;
        self->prev_lp_cutoff = current_lp_cutoff;
        self->prev_lp_resonance = current_lp_resonance;
        self->prev_hp_slope = (float)hp_slope_selector;
        self->prev_lp_slope = (float)lp_slope_selector;
    }

    // Process audio through filters - BOTH CHANNELS TOGETHER in same loop
    for (uint32_t i = 0; i < n_samples; i++) {
        // Load both channels
        float left_signal = self->input_left[i];
        float right_signal = self->input_right[i];
        
        // Process high-pass filter on both channels (skip if HPF is off)
        if (hp_slope_selector != 5) {
            switch (hp_slope_selector) {
                case 4:  // -48dB Butterworth HPF
                    for (int j = 0; j < 4; ++j) {
                        // STATE OPTIMIZATION: Use temporary variables
                        float x1_left = self->hp48_bw_x1[0][j];
                        float x2_left = self->hp48_bw_x2[0][j];
                        float y1_left = self->hp48_bw_y1[0][j];
                        float y2_left = self->hp48_bw_y2[0][j];
                        
                        float x1_right = self->hp48_bw_x1[1][j];
                        float x2_right = self->hp48_bw_x2[1][j];
                        float y1_right = self->hp48_bw_y1[1][j];
                        float y2_right = self->hp48_bw_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->hp48_bw_b0[j] * left_signal + self->hp48_bw_b1[j] * x1_left + 
                                            self->hp48_bw_b2[j] * x2_left - self->hp48_bw_a1[j] * y1_left - 
                                            self->hp48_bw_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->hp48_bw_b0[j] * right_signal + self->hp48_bw_b1[j] * x1_right + 
                                             self->hp48_bw_b2[j] * x2_right - self->hp48_bw_a1[j] * y1_right - 
                                             self->hp48_bw_a2[j] * y2_right;
                        
                        // Update state variables (minimal memory writes)
                        self->hp48_bw_x2[0][j] = x1_left;
                        self->hp48_bw_x1[0][j] = left_signal;
                        self->hp48_bw_y2[0][j] = y1_left;
                        self->hp48_bw_y1[0][j] = temp_out_left;
                        
                        self->hp48_bw_x2[1][j] = x1_right;
                        self->hp48_bw_x1[1][j] = right_signal;
                        self->hp48_bw_y2[1][j] = y1_right;
                        self->hp48_bw_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 3:  // -24dB Butterworth HPF
                    for (int j = 0; j < 2; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->hp24_bw_x1[0][j];
                        float x2_left = self->hp24_bw_x2[0][j];
                        float y1_left = self->hp24_bw_y1[0][j];
                        float y2_left = self->hp24_bw_y2[0][j];
                        
                        float x1_right = self->hp24_bw_x1[1][j];
                        float x2_right = self->hp24_bw_x2[1][j];
                        float y1_right = self->hp24_bw_y1[1][j];
                        float y2_right = self->hp24_bw_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->hp24_bw_b0[j] * left_signal + self->hp24_bw_b1[j] * x1_left + 
                                            self->hp24_bw_b2[j] * x2_left - self->hp24_bw_a1[j] * y1_left - 
                                            self->hp24_bw_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->hp24_bw_b0[j] * right_signal + self->hp24_bw_b1[j] * x1_right + 
                                             self->hp24_bw_b2[j] * x2_right - self->hp24_bw_a1[j] * y1_right - 
                                             self->hp24_bw_a2[j] * y2_right;
                        
                        // Update state variables
                        self->hp24_bw_x2[0][j] = x1_left;
                        self->hp24_bw_x1[0][j] = left_signal;
                        self->hp24_bw_y2[0][j] = y1_left;
                        self->hp24_bw_y1[0][j] = temp_out_left;
                        
                        self->hp24_bw_x2[1][j] = x1_right;
                        self->hp24_bw_x1[1][j] = right_signal;
                        self->hp24_bw_y2[1][j] = y1_right;
                        self->hp24_bw_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 2:  // -48dB variable HPF
                    for (int j = 0; j < 4; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->hp48_x1[0][j];
                        float x2_left = self->hp48_x2[0][j];
                        float y1_left = self->hp48_y1[0][j];
                        float y2_left = self->hp48_y2[0][j];
                        
                        float x1_right = self->hp48_x1[1][j];
                        float x2_right = self->hp48_x2[1][j];
                        float y1_right = self->hp48_y1[1][j];
                        float y2_right = self->hp48_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->hp48_b0[j] * left_signal + self->hp48_b1[j] * x1_left + 
                                            self->hp48_b2[j] * x2_left - self->hp48_a1[j] * y1_left - 
                                            self->hp48_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->hp48_b0[j] * right_signal + self->hp48_b1[j] * x1_right + 
                                             self->hp48_b2[j] * x2_right - self->hp48_a1[j] * y1_right - 
                                             self->hp48_a2[j] * y2_right;
                        
                        // Update state variables
                        self->hp48_x2[0][j] = x1_left;
                        self->hp48_x1[0][j] = left_signal;
                        self->hp48_y2[0][j] = y1_left;
                        self->hp48_y1[0][j] = temp_out_left;
                        
                        self->hp48_x2[1][j] = x1_right;
                        self->hp48_x1[1][j] = right_signal;
                        self->hp48_y2[1][j] = y1_right;
                        self->hp48_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 1:  // -24dB variable HPF
                    for (int j = 0; j < 2; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->hp24_x1[0][j];
                        float x2_left = self->hp24_x2[0][j];
                        float y1_left = self->hp24_y1[0][j];
                        float y2_left = self->hp24_y2[0][j];
                        
                        float x1_right = self->hp24_x1[1][j];
                        float x2_right = self->hp24_x2[1][j];
                        float y1_right = self->hp24_y1[1][j];
                        float y2_right = self->hp24_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->hp24_b0[j] * left_signal + self->hp24_b1[j] * x1_left + 
                                            self->hp24_b2[j] * x2_left - self->hp24_a1[j] * y1_left - 
                                            self->hp24_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->hp24_b0[j] * right_signal + self->hp24_b1[j] * x1_right + 
                                             self->hp24_b2[j] * x2_right - self->hp24_a1[j] * y1_right - 
                                             self->hp24_a2[j] * y2_right;
                        
                        // Update state variables
                        self->hp24_x2[0][j] = x1_left;
                        self->hp24_x1[0][j] = left_signal;
                        self->hp24_y2[0][j] = y1_left;
                        self->hp24_y1[0][j] = temp_out_left;
                        
                        self->hp24_x2[1][j] = x1_right;
                        self->hp24_x1[1][j] = right_signal;
                        self->hp24_y2[1][j] = y1_right;
                        self->hp24_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 0:  // -12dB variable HPF
                default:
                    // STATE OPTIMIZATION
                    float x1_left = self->hp12_x1[0];
                    float x2_left = self->hp12_x2[0];
                    float y1_left = self->hp12_y1[0];
                    float y2_left = self->hp12_y2[0];
                    
                    float x1_right = self->hp12_x1[1];
                    float x2_right = self->hp12_x2[1];
                    float y1_right = self->hp12_y1[1];
                    float y2_right = self->hp12_y2[1];
                    
                    // Process left channel
                    float temp_out_left = self->hp12_b0 * left_signal + self->hp12_b1 * x1_left + 
                                        self->hp12_b2 * x2_left - self->hp12_a1 * y1_left - 
                                        self->hp12_a2 * y2_left;
                    
                    // Process right channel
                    float temp_out_right = self->hp12_b0 * right_signal + self->hp12_b1 * x1_right + 
                                         self->hp12_b2 * x2_right - self->hp12_a1 * y1_right - 
                                         self->hp12_a2 * y2_right;
                    
                    // Update state variables
                    self->hp12_x2[0] = x1_left;
                    self->hp12_x1[0] = left_signal;
                    self->hp12_y2[0] = y1_left;
                    self->hp12_y1[0] = temp_out_left;
                    
                    self->hp12_x2[1] = x1_right;
                    self->hp12_x1[1] = right_signal;
                    self->hp12_y2[1] = y1_right;
                    self->hp12_y1[1] = temp_out_right;
                    
                    left_signal = temp_out_left;
                    right_signal = temp_out_right;
                    break;
            }
        }
        
        // Process low-pass filter on both channels (skip if LPF is off)
        if (lp_slope_selector != 5) {
            switch (lp_slope_selector) {
                case 4:  // -48dB Butterworth LPF
                    for (int j = 0; j < 4; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->lp48_bw_x1[0][j];
                        float x2_left = self->lp48_bw_x2[0][j];
                        float y1_left = self->lp48_bw_y1[0][j];
                        float y2_left = self->lp48_bw_y2[0][j];
                        
                        float x1_right = self->lp48_bw_x1[1][j];
                        float x2_right = self->lp48_bw_x2[1][j];
                        float y1_right = self->lp48_bw_y1[1][j];
                        float y2_right = self->lp48_bw_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->lp48_bw_b0[j] * left_signal + self->lp48_bw_b1[j] * x1_left + 
                                            self->lp48_bw_b2[j] * x2_left - self->lp48_bw_a1[j] * y1_left - 
                                            self->lp48_bw_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->lp48_bw_b0[j] * right_signal + self->lp48_bw_b1[j] * x1_right + 
                                             self->lp48_bw_b2[j] * x2_right - self->lp48_bw_a1[j] * y1_right - 
                                             self->lp48_bw_a2[j] * y2_right;
                        
                        // Update state variables
                        self->lp48_bw_x2[0][j] = x1_left;
                        self->lp48_bw_x1[0][j] = left_signal;
                        self->lp48_bw_y2[0][j] = y1_left;
                        self->lp48_bw_y1[0][j] = temp_out_left;
                        
                        self->lp48_bw_x2[1][j] = x1_right;
                        self->lp48_bw_x1[1][j] = right_signal;
                        self->lp48_bw_y2[1][j] = y1_right;
                        self->lp48_bw_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 3:  // -24dB Butterworth LPF
                    for (int j = 0; j < 2; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->lp24_bw_x1[0][j];
                        float x2_left = self->lp24_bw_x2[0][j];
                        float y1_left = self->lp24_bw_y1[0][j];
                        float y2_left = self->lp24_bw_y2[0][j];
                        
                        float x1_right = self->lp24_bw_x1[1][j];
                        float x2_right = self->lp24_bw_x2[1][j];
                        float y1_right = self->lp24_bw_y1[1][j];
                        float y2_right = self->lp24_bw_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->lp24_bw_b0[j] * left_signal + self->lp24_bw_b1[j] * x1_left + 
                                            self->lp24_bw_b2[j] * x2_left - self->lp24_bw_a1[j] * y1_left - 
                                            self->lp24_bw_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->lp24_bw_b0[j] * right_signal + self->lp24_bw_b1[j] * x1_right + 
                                             self->lp24_bw_b2[j] * x2_right - self->lp24_bw_a1[j] * y1_right - 
                                             self->lp24_bw_a2[j] * y2_right;
                        
                        // Update state variables
                        self->lp24_bw_x2[0][j] = x1_left;
                        self->lp24_bw_x1[0][j] = left_signal;
                        self->lp24_bw_y2[0][j] = y1_left;
                        self->lp24_bw_y1[0][j] = temp_out_left;
                        
                        self->lp24_bw_x2[1][j] = x1_right;
                        self->lp24_bw_x1[1][j] = right_signal;
                        self->lp24_bw_y2[1][j] = y1_right;
                        self->lp24_bw_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 2:  // -48dB variable LPF
                    for (int j = 0; j < 4; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->lp48_x1[0][j];
                        float x2_left = self->lp48_x2[0][j];
                        float y1_left = self->lp48_y1[0][j];
                        float y2_left = self->lp48_y2[0][j];
                        
                        float x1_right = self->lp48_x1[1][j];
                        float x2_right = self->lp48_x2[1][j];
                        float y1_right = self->lp48_y1[1][j];
                        float y2_right = self->lp48_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->lp48_b0[j] * left_signal + self->lp48_b1[j] * x1_left + 
                                            self->lp48_b2[j] * x2_left - self->lp48_a1[j] * y1_left - 
                                            self->lp48_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->lp48_b0[j] * right_signal + self->lp48_b1[j] * x1_right + 
                                             self->lp48_b2[j] * x2_right - self->lp48_a1[j] * y1_right - 
                                             self->lp48_a2[j] * y2_right;
                        
                        // Update state variables
                        self->lp48_x2[0][j] = x1_left;
                        self->lp48_x1[0][j] = left_signal;
                        self->lp48_y2[0][j] = y1_left;
                        self->lp48_y1[0][j] = temp_out_left;
                        
                        self->lp48_x2[1][j] = x1_right;
                        self->lp48_x1[1][j] = right_signal;
                        self->lp48_y2[1][j] = y1_right;
                        self->lp48_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 1:  // -24dB variable LPF
                    for (int j = 0; j < 2; ++j) {
                        // STATE OPTIMIZATION
                        float x1_left = self->lp24_x1[0][j];
                        float x2_left = self->lp24_x2[0][j];
                        float y1_left = self->lp24_y1[0][j];
                        float y2_left = self->lp24_y2[0][j];
                        
                        float x1_right = self->lp24_x1[1][j];
                        float x2_right = self->lp24_x2[1][j];
                        float y1_right = self->lp24_y1[1][j];
                        float y2_right = self->lp24_y2[1][j];
                        
                        // Process left channel
                        float temp_out_left = self->lp24_b0[j] * left_signal + self->lp24_b1[j] * x1_left + 
                                            self->lp24_b2[j] * x2_left - self->lp24_a1[j] * y1_left - 
                                            self->lp24_a2[j] * y2_left;
                        
                        // Process right channel
                        float temp_out_right = self->lp24_b0[j] * right_signal + self->lp24_b1[j] * x1_right + 
                                             self->lp24_b2[j] * x2_right - self->lp24_a1[j] * y1_right - 
                                             self->lp24_a2[j] * y2_right;
                        
                        // Update state variables
                        self->lp24_x2[0][j] = x1_left;
                        self->lp24_x1[0][j] = left_signal;
                        self->lp24_y2[0][j] = y1_left;
                        self->lp24_y1[0][j] = temp_out_left;
                        
                        self->lp24_x2[1][j] = x1_right;
                        self->lp24_x1[1][j] = right_signal;
                        self->lp24_y2[1][j] = y1_right;
                        self->lp24_y1[1][j] = temp_out_right;
                        
                        left_signal = temp_out_left;
                        right_signal = temp_out_right;
                    }
                    break;
                    
                case 0:  // -12dB variable LPF
                default:
                    // STATE OPTIMIZATION
                    float x1_left = self->lp12_x1[0];
                    float x2_left = self->lp12_x2[0];
                    float y1_left = self->lp12_y1[0];
                    float y2_left = self->lp12_y2[0];
                    
                    float x1_right = self->lp12_x1[1];
                    float x2_right = self->lp12_x2[1];
                    float y1_right = self->lp12_y1[1];
                    float y2_right = self->lp12_y2[1];
                    
                    // Process left channel
                    float temp_out_left = self->lp12_b0 * left_signal + self->lp12_b1 * x1_left + 
                                        self->lp12_b2 * x2_left - self->lp12_a1 * y1_left - 
                                        self->lp12_a2 * y2_left;
                    
                    // Process right channel
                    float temp_out_right = self->lp12_b0 * right_signal + self->lp12_b1 * x1_right + 
                                         self->lp12_b2 * x2_right - self->lp12_a1 * y1_right - 
                                         self->lp12_a2 * y2_right;
                    
                    // Update state variables
                    self->lp12_x2[0] = x1_left;
                    self->lp12_x1[0] = left_signal;
                    self->lp12_y2[0] = y1_left;
                    self->lp12_y1[0] = temp_out_left;
                    
                    self->lp12_x2[1] = x1_right;
                    self->lp12_x1[1] = right_signal;
                    self->lp12_y2[1] = y1_right;
                    self->lp12_y1[1] = temp_out_right;
                    
                    left_signal = temp_out_left;
                    right_signal = temp_out_right;
                    break;
            }
        }
        // Store both channels
        self->output_left[i] = left_signal;
        self->output_right[i] = right_signal;
    }
}

static void cleanup(LV2_Handle instance) {
    free(instance);
}

static const LV2_Descriptor descriptor = {
    DILOPHILTER_URI,
    instantiate,
    connect_port,
    NULL,
    run,
    NULL,
    cleanup,
    NULL
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index) {
    return (index == 0) ? &descriptor : NULL;
}
