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
#include <string.h>

#define DILOPHILTER_URI "urn:simdott:dilophilter"

typedef enum {
    INPUT_LEFT = 0,
    INPUT_RIGHT = 1,
    OUTPUT_LEFT = 2,
    OUTPUT_RIGHT = 3,
    HP_SLOPE = 4,
    HP_CUTOFF = 5,
    HP_RESONANCE = 6,
    HP_GAIN = 7,
    LP_SLOPE = 8,
    LP_CUTOFF = 9,
    LP_RESONANCE = 10,
    LP_GAIN = 11,
} PortIndex;

typedef struct {
    const float* input_left;
    const float* input_right;
    const float* hp_slope;
    const float* hp_cutoff;
    const float* hp_resonance;
    const float* hp_gain;
    const float* lp_slope;
    const float* lp_cutoff;
    const float* lp_resonance;
    const float* lp_gain;
    float* output_left;
    float* output_right;
    
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
    
    float hp48_lr_a1[4], hp48_lr_a2[4], hp48_lr_b0[4], hp48_lr_b1[4], hp48_lr_b2[4];
    float hp48_lr_x1[2][4], hp48_lr_x2[2][4], hp48_lr_y1[2][4], hp48_lr_y2[2][4];
    
    float hp96_lr_a1[8], hp96_lr_a2[8], hp96_lr_b0[8], hp96_lr_b1[8], hp96_lr_b2[8];
    float hp96_lr_x1[2][8], hp96_lr_x2[2][8], hp96_lr_y1[2][8], hp96_lr_y2[2][8];
    
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
    
    float lp48_lr_a1[4], lp48_lr_a2[4], lp48_lr_b0[4], lp48_lr_b1[4], lp48_lr_b2[4];
    float lp48_lr_x1[2][4], lp48_lr_x2[2][4], lp48_lr_y1[2][4], lp48_lr_y2[2][4];
    
    float lp96_lr_a1[8], lp96_lr_a2[8], lp96_lr_b0[8], lp96_lr_b1[8], lp96_lr_b2[8];
    float lp96_lr_x1[2][8], lp96_lr_x2[2][8], lp96_lr_y1[2][8], lp96_lr_y2[2][8];

    float ls_a1, ls_a2, ls_b0, ls_b1, ls_b2;
    float ls_x1[2], ls_x2[2], ls_y1[2], ls_y2[2];
    
    float hs_a1, hs_a2, hs_b0, hs_b1, hs_b2;
    float hs_x1[2], hs_x2[2], hs_y1[2], hs_y2[2];

    float prev_hp_cutoff, prev_hp_resonance, prev_hp_gain;
    float prev_lp_cutoff, prev_lp_resonance, prev_lp_gain;
    float sample_rate;
    float prev_hp_slope, prev_lp_slope;

} Dilophilter;

static float linear_gain(float db) {
    if (db <= -120.0f) return 0.0f;
    return powf(10.0f, db / 20.0f);
}

// Convert relative resonance (0-100%) to Q value based on filter type
static float relative_resonance_to_q_12db(float relative_percent) {
    // Convert from 0-100 to 0-1 range
    float relative = relative_percent / 100.0f;
    
    // -12dB: 0% = Q=0.126, 50% = Q=0.707, 100% = Q=4.0
    if (relative <= 0.5f) {
        float t = relative / 0.5f;
        return 0.126f + (0.707f - 0.126f) * t;
    } else {
        float t = (relative - 0.5f) / 0.5f;
        return 0.707f + (4.0f - 0.707f) * t;
    }
}

static float relative_resonance_to_q_24db(float relative_percent) {
    float relative = relative_percent / 100.0f;
    
    // -24dB: 0% = Q=0.356, 50% = Q=0.707, 100% = Q=2.0
    if (relative <= 0.5f) {
        float t = relative / 0.5f;
        return 0.356f + (0.707f - 0.356f) * t;
    } else {
        float t = (relative - 0.5f) / 0.5f;
        return 0.707f + (2.0f - 0.707f) * t;
    }
}

static float relative_resonance_to_q_48db_classic(float relative_percent) {
    float relative = relative_percent / 100.0f;
    float Q;
    
    // Classic 48dB: 0% = Q=0.597, 50% = Q=0.7071, 100% = Q=1.414
    if (relative <= 0.5f) {
        float t = relative / 0.5f;
        Q = 0.597f + (0.7071f - 0.597f) * t;
    } else {
        float t = (relative - 0.5f) / 0.5f;
        Q = 0.7071f + (1.414f - 0.7071f) * t;
    }
    
    return Q;
}

static float relative_resonance_to_q_shelf(float relative_percent) {
    float relative = relative_percent / 100.0f;
    
    // Shelf filters: 0% = Q=0.3, 50% = Q=0.707, 100% = Q=1.5
    if (relative <= 0.5f) {
        float t = relative / 0.5f;
        return 0.3f + (0.707f - 0.3f) * t;
    } else {
        float t = (relative - 0.5f) / 0.5f;
        return 0.707f + (1.5f - 0.707f) * t;
    }
}

static void calculate_low_shelf_coefficients(Dilophilter* self, float cutoff, float relative_resonance, float gain_db) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;
    
    float Q = relative_resonance_to_q_shelf(relative_resonance);
    if (Q < 0.3f) Q = 0.3f;
    if (Q > 1.5f) Q = 1.5f;
    
    float gain_linear = linear_gain(gain_db);
    float A = sqrtf(gain_linear);
    
    float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    float sin_omega = sinf(omega);
    float cos_omega = cosf(omega);
    float alpha = sin_omega / (2.0f * Q);
    float sqrtA = sqrtf(A);
    
    float b0 = A * ((A + 1.0f) - (A - 1.0f) * cos_omega + 2.0f * sqrtA * alpha);
    float b1 = 2.0f * A * ((A - 1.0f) - (A + 1.0f) * cos_omega);
    float b2 = A * ((A + 1.0f) - (A - 1.0f) * cos_omega - 2.0f * sqrtA * alpha);
    float a0 = (A + 1.0f) + (A - 1.0f) * cos_omega + 2.0f * sqrtA * alpha;
    float a1 = -2.0f * ((A - 1.0f) + (A + 1.0f) * cos_omega);
    float a2 = (A + 1.0f) + (A - 1.0f) * cos_omega - 2.0f * sqrtA * alpha;
    
    self->ls_b0 = b0 / a0;
    self->ls_b1 = b1 / a0;
    self->ls_b2 = b2 / a0;
    self->ls_a1 = a1 / a0;
    self->ls_a2 = a2 / a0;
}

static void calculate_high_shelf_coefficients(Dilophilter* self, float cutoff, float relative_resonance, float gain_db) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;
    
    float Q = relative_resonance_to_q_shelf(relative_resonance);
    if (Q < 0.3f) Q = 0.3f;
    if (Q > 1.5f) Q = 1.5f;
    
    float gain_linear = linear_gain(gain_db);
    float A = sqrtf(gain_linear);
    
    float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    float sin_omega = sinf(omega);
    float cos_omega = cosf(omega);
    float alpha = sin_omega / (2.0f * Q);
    float sqrtA = sqrtf(A);
    
    float b0 = A * ((A + 1.0f) + (A - 1.0f) * cos_omega + 2.0f * sqrtA * alpha);
    float b1 = -2.0f * A * ((A - 1.0f) + (A + 1.0f) * cos_omega);
    float b2 = A * ((A + 1.0f) + (A - 1.0f) * cos_omega - 2.0f * sqrtA * alpha);
    float a0 = (A + 1.0f) - (A - 1.0f) * cos_omega + 2.0f * sqrtA * alpha;
    float a1 = 2.0f * ((A - 1.0f) - (A + 1.0f) * cos_omega);
    float a2 = (A + 1.0f) - (A - 1.0f) * cos_omega - 2.0f * sqrtA * alpha;
    
    self->hs_b0 = b0 / a0;
    self->hs_b1 = b1 / a0;
    self->hs_b2 = b2 / a0;
    self->hs_a1 = a1 / a0;
    self->hs_a2 = a2 / a0;
}

static void calculate_hp12_coefficients(Dilophilter* self, float cutoff, float relative_resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;
    
    float Q = relative_resonance_to_q_12db(relative_resonance);
    if (Q < 0.126f) Q = 0.126f;
    if (Q > 4.0f) Q = 4.0f;

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

static void calculate_hp24_coefficients(Dilophilter* self, float cutoff, float relative_resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    float Q = relative_resonance_to_q_24db(relative_resonance);
    if (Q < 0.356f) Q = 0.356f;
    if (Q > 2.0f) Q = 2.0f;

    for (int i = 0; i < 2; ++i) {
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

static void calculate_hp48_coefficients(Dilophilter* self, float cutoff, float relative_resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    float Q = relative_resonance_to_q_48db_classic(relative_resonance);
    if (Q < 0.597f) Q = 0.597f;
    if (Q > 1.414f) Q = 1.414f;

    for (int i = 0; i < 4; ++i) {
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

static void calculate_hp48_lr_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    const float Q_values[4] = {0.5412f, 1.3066f, 0.5412f, 1.3066f};

    for (int i = 0; i < 4; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        const float b0 = (1.0f + cos_omega) / 2.0f;
        const float b1 = -(1.0f + cos_omega);
        const float b2 = (1.0f + cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->hp48_lr_b0[i] = b0 / a0;
        self->hp48_lr_b1[i] = b1 / a0;
        self->hp48_lr_b2[i] = b2 / a0;
        self->hp48_lr_a1[i] = a1 / a0;
        self->hp48_lr_a2[i] = a2 / a0;
    }
}

static void calculate_hp96_lr_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    const float Q_values[8] = {0.5098f, 0.6013f, 0.8999f, 2.5628f, 0.5098f, 0.6013f, 0.8999f, 2.5628f};

    for (int i = 0; i < 8; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        const float b0 = (1.0f + cos_omega) / 2.0f;
        const float b1 = -(1.0f + cos_omega);
        const float b2 = (1.0f + cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->hp96_lr_b0[i] = b0 / a0;
        self->hp96_lr_b1[i] = b1 / a0;
        self->hp96_lr_b2[i] = b2 / a0;
        self->hp96_lr_a1[i] = a1 / a0;
        self->hp96_lr_a2[i] = a2 / a0;
    }
}

static void calculate_lp12_coefficients(Dilophilter* self, float cutoff, float relative_resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;
    
    float Q = relative_resonance_to_q_12db(relative_resonance);
    if (Q < 0.126f) Q = 0.126f;
    if (Q > 4.0f) Q = 4.0f;

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

static void calculate_lp24_coefficients(Dilophilter* self, float cutoff, float relative_resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    float Q = relative_resonance_to_q_24db(relative_resonance);
    if (Q < 0.356f) Q = 0.356f;
    if (Q > 2.0f) Q = 2.0f;

    for (int i = 0; i < 2; ++i) {
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

static void calculate_lp48_coefficients(Dilophilter* self, float cutoff, float relative_resonance) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    
    float Q = relative_resonance_to_q_48db_classic(relative_resonance);
    if (Q < 0.597f) Q = 0.597f;
    if (Q > 1.414f) Q = 1.414f;

    for (int i = 0; i < 4; ++i) {
        const float alpha = sin_omega / (2.0f * Q);
        const float a0 = 1.0f + alpha;
        // Low-pass coefficients (note: b0, b1, b2 are different from HP)
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

static void calculate_lp48_lr_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    const float Q_values[4] = {0.5412f, 1.3066f, 0.5412f, 1.3066f};

    for (int i = 0; i < 4; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        const float b0 = (1.0f - cos_omega) / 2.0f;
        const float b1 = 1.0f - cos_omega;
        const float b2 = (1.0f - cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->lp48_lr_b0[i] = b0 / a0;
        self->lp48_lr_b1[i] = b1 / a0;
        self->lp48_lr_b2[i] = b2 / a0;
        self->lp48_lr_a1[i] = a1 / a0;
        self->lp48_lr_a2[i] = a2 / a0;
    }
}

static void calculate_lp96_lr_coefficients(Dilophilter* self, float cutoff) {
    if (cutoff < 20.0f) cutoff = 20.0f;
    if (cutoff > 20000.0f) cutoff = 20000.0f;

    const float omega = 2.0f * M_PI * cutoff / self->sample_rate;
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    const float Q_values[8] = {0.5098f, 0.6013f, 0.8999f, 2.5628f, 0.5098f, 0.6013f, 0.8999f, 2.5628f};

    for (int i = 0; i < 8; ++i) {
        const float alpha = sin_omega / (2.0f * Q_values[i]);
        const float a0 = 1.0f + alpha;
        const float b0 = (1.0f - cos_omega) / 2.0f;
        const float b1 = 1.0f - cos_omega;
        const float b2 = (1.0f - cos_omega) / 2.0f;
        const float a1 = -2.0f * cos_omega;
        const float a2 = 1.0f - alpha;

        self->lp96_lr_b0[i] = b0 / a0;
        self->lp96_lr_b1[i] = b1 / a0;
        self->lp96_lr_b2[i] = b2 / a0;
        self->lp96_lr_a1[i] = a1 / a0;
        self->lp96_lr_a2[i] = a2 / a0;
    }
}

static LV2_Handle instantiate(
    const LV2_Descriptor* descriptor,
    double sample_rate,
    const char* bundle_path,
    const LV2_Feature* const* features) {
    
    Dilophilter* self = (Dilophilter*)calloc(1, sizeof(Dilophilter));
    if (!self) return NULL;
    
    self->sample_rate = (float)sample_rate;
    self->prev_hp_cutoff = 30.0f;
    self->prev_hp_resonance = 50.0f;  // 50%
    self->prev_hp_gain = 0.0f;
    self->prev_lp_cutoff = 10000.0f;
    self->prev_lp_resonance = 50.0f;  // 50%
    self->prev_lp_gain = 0.0f;
    self->prev_hp_slope = 0.0f;
    self->prev_lp_slope = 0.0f;
    
    calculate_low_shelf_coefficients(self, 30.0f, 50.0f, 0.0f);
    calculate_high_shelf_coefficients(self, 10000.0f, 50.0f, 0.0f);
    
    calculate_hp12_coefficients(self, 30.0f, 50.0f);
    calculate_hp24_coefficients(self, 30.0f, 50.0f);
    calculate_hp48_coefficients(self, 30.0f, 50.0f);
    calculate_hp24_bw_coefficients(self, 30.0f);
    calculate_hp48_bw_coefficients(self, 30.0f);
    calculate_hp48_lr_coefficients(self, 30.0f);
    calculate_hp96_lr_coefficients(self, 30.0f);

    calculate_lp12_coefficients(self, 10000.0f, 50.0f);
    calculate_lp24_coefficients(self, 10000.0f, 50.0f);
    calculate_lp48_coefficients(self, 10000.0f, 50.0f);
    calculate_lp24_bw_coefficients(self, 10000.0f);
    calculate_lp48_bw_coefficients(self, 10000.0f);
    calculate_lp48_lr_coefficients(self, 10000.0f);
    calculate_lp96_lr_coefficients(self, 10000.0f);
    
    return (LV2_Handle)self;
}

static void connect_port(LV2_Handle instance, uint32_t port, void* data) {
    Dilophilter* self = (Dilophilter*)instance;
    
    switch ((PortIndex)port) {
        case INPUT_LEFT: self->input_left = (const float*)data; break;
        case INPUT_RIGHT: self->input_right = (const float*)data; break;
        case OUTPUT_LEFT: self->output_left = (float*)data; break;
        case OUTPUT_RIGHT: self->output_right = (float*)data; break;
        case HP_SLOPE: self->hp_slope = (const float*)data; break;
        case HP_CUTOFF: self->hp_cutoff = (const float*)data; break;
        case HP_RESONANCE: self->hp_resonance = (const float*)data; break;
        case HP_GAIN: self->hp_gain = (const float*)data; break;
        case LP_SLOPE: self->lp_slope = (const float*)data; break;
        case LP_CUTOFF: self->lp_cutoff = (const float*)data; break;
        case LP_RESONANCE: self->lp_resonance = (const float*)data; break;
        case LP_GAIN: self->lp_gain = (const float*)data; break;
    }
}

static void run(LV2_Handle instance, uint32_t n_samples) {
    Dilophilter* self = (Dilophilter*)instance;
    
    int hp_slope_selector = (int)*(self->hp_slope);
    int lp_slope_selector = (int)*(self->lp_slope);
    float current_hp_cutoff = *(self->hp_cutoff);
    float current_hp_resonance = *(self->hp_resonance);
    float current_hp_gain = *(self->hp_gain);
    float current_lp_cutoff = *(self->lp_cutoff);
    float current_lp_resonance = *(self->lp_resonance);
    float current_lp_gain = *(self->lp_gain);
    
    if (hp_slope_selector == 0 && lp_slope_selector == 0) {
        for (uint32_t i = 0; i < n_samples; i++) {
            self->output_left[i] = self->input_left[i];
            self->output_right[i] = self->input_right[i];
        }
        return;
    }
    
    int parameters_changed = (current_hp_cutoff != self->prev_hp_cutoff) ||
                         (current_lp_cutoff != self->prev_lp_cutoff) ||
                         (current_hp_resonance != self->prev_hp_resonance) ||
                         (current_lp_resonance != self->prev_lp_resonance) ||
                         (hp_slope_selector != (int)self->prev_hp_slope) ||
                         (lp_slope_selector != (int)self->prev_lp_slope);

    if ((hp_slope_selector == 7) || (hp_slope_selector == 8)) {
        parameters_changed |= (current_hp_gain != self->prev_hp_gain);
    }
    if ((lp_slope_selector == 7) || (lp_slope_selector == 8)) {
        parameters_changed |= (current_lp_gain != self->prev_lp_gain);
    }

    if (parameters_changed) {
        if (hp_slope_selector != 0) {
            if (hp_slope_selector == 8) {
                calculate_low_shelf_coefficients(self, current_hp_cutoff, current_hp_resonance, current_hp_gain);
            } else {
                if (hp_slope_selector == 1) calculate_hp12_coefficients(self, current_hp_cutoff, current_hp_resonance);
                else if (hp_slope_selector == 2) calculate_hp24_coefficients(self, current_hp_cutoff, current_hp_resonance);
                else if (hp_slope_selector == 3) calculate_hp48_coefficients(self, current_hp_cutoff, current_hp_resonance);
                else if (hp_slope_selector == 4) calculate_hp24_bw_coefficients(self, current_hp_cutoff);
                else if (hp_slope_selector == 5) calculate_hp48_bw_coefficients(self, current_hp_cutoff);
                else if (hp_slope_selector == 6) calculate_hp48_lr_coefficients(self, current_hp_cutoff);
                else if (hp_slope_selector == 7) calculate_hp96_lr_coefficients(self, current_hp_cutoff);
            }
        }
        
        if (lp_slope_selector != 0) {
            if (lp_slope_selector == 8) {
                calculate_high_shelf_coefficients(self, current_lp_cutoff, current_lp_resonance, current_lp_gain);
            } else {
                if (lp_slope_selector == 1) calculate_lp12_coefficients(self, current_lp_cutoff, current_lp_resonance);
                else if (lp_slope_selector == 2) calculate_lp24_coefficients(self, current_lp_cutoff, current_lp_resonance);
                else if (lp_slope_selector == 3) calculate_lp48_coefficients(self, current_lp_cutoff, current_lp_resonance);
                else if (lp_slope_selector == 4) calculate_lp24_bw_coefficients(self, current_lp_cutoff);
                else if (lp_slope_selector == 5) calculate_lp48_bw_coefficients(self, current_lp_cutoff);
                else if (lp_slope_selector == 6) calculate_lp48_lr_coefficients(self, current_lp_cutoff);
                else if (lp_slope_selector == 7) calculate_lp96_lr_coefficients(self, current_lp_cutoff);
            }
        }
        
        self->prev_hp_cutoff = current_hp_cutoff;
        self->prev_lp_cutoff = current_lp_cutoff;
        self->prev_hp_resonance = current_hp_resonance;
        self->prev_lp_resonance = current_lp_resonance;
        self->prev_hp_slope = (float)hp_slope_selector;
        self->prev_lp_slope = (float)lp_slope_selector;
        
        if ((hp_slope_selector == 7) || (hp_slope_selector == 8)) {
            self->prev_hp_gain = current_hp_gain;
        }
        if ((lp_slope_selector == 7) || (lp_slope_selector == 8)) {
            self->prev_lp_gain = current_lp_gain;
        }
    }

    for (uint32_t i = 0; i < n_samples; i++) {
        float left_signal = self->input_left[i];
        float right_signal = self->input_right[i];
        
        if (hp_slope_selector != 0) {
            if (hp_slope_selector == 8) {
                float x1L = self->ls_x1[0], x2L = self->ls_x2[0], y1L = self->ls_y1[0], y2L = self->ls_y2[0];
                float x1R = self->ls_x1[1], x2R = self->ls_x2[1], y1R = self->ls_y1[1], y2R = self->ls_y2[1];
                
                float outL = self->ls_b0 * left_signal + self->ls_b1 * x1L + self->ls_b2 * x2L - self->ls_a1 * y1L - self->ls_a2 * y2L;
                float outR = self->ls_b0 * right_signal + self->ls_b1 * x1R + self->ls_b2 * x2R - self->ls_a1 * y1R - self->ls_a2 * y2R;
                
                self->ls_x2[0] = x1L; self->ls_x1[0] = left_signal; self->ls_y2[0] = y1L; self->ls_y1[0] = outL;
                self->ls_x2[1] = x1R; self->ls_x1[1] = right_signal; self->ls_y2[1] = y1R; self->ls_y1[1] = outR;
                
                left_signal = outL; right_signal = outR;
            } else if (hp_slope_selector == 7) {
                for (int j = 0; j < 8; ++j) {
                    float x1L = self->hp96_lr_x1[0][j], x2L = self->hp96_lr_x2[0][j], y1L = self->hp96_lr_y1[0][j], y2L = self->hp96_lr_y2[0][j];
                    float x1R = self->hp96_lr_x1[1][j], x2R = self->hp96_lr_x2[1][j], y1R = self->hp96_lr_y1[1][j], y2R = self->hp96_lr_y2[1][j];
                    
                    float outL = self->hp96_lr_b0[j] * left_signal + self->hp96_lr_b1[j] * x1L + self->hp96_lr_b2[j] * x2L - self->hp96_lr_a1[j] * y1L - self->hp96_lr_a2[j] * y2L;
                    float outR = self->hp96_lr_b0[j] * right_signal + self->hp96_lr_b1[j] * x1R + self->hp96_lr_b2[j] * x2R - self->hp96_lr_a1[j] * y1R - self->hp96_lr_a2[j] * y2R;
                    
                    self->hp96_lr_x2[0][j] = x1L; self->hp96_lr_x1[0][j] = left_signal; self->hp96_lr_y2[0][j] = y1L; self->hp96_lr_y1[0][j] = outL;
                    self->hp96_lr_x2[1][j] = x1R; self->hp96_lr_x1[1][j] = right_signal; self->hp96_lr_y2[1][j] = y1R; self->hp96_lr_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (hp_slope_selector == 6) {
                for (int j = 0; j < 4; ++j) {
                    float x1L = self->hp48_lr_x1[0][j], x2L = self->hp48_lr_x2[0][j], y1L = self->hp48_lr_y1[0][j], y2L = self->hp48_lr_y2[0][j];
                    float x1R = self->hp48_lr_x1[1][j], x2R = self->hp48_lr_x2[1][j], y1R = self->hp48_lr_y1[1][j], y2R = self->hp48_lr_y2[1][j];
                    
                    float outL = self->hp48_lr_b0[j] * left_signal + self->hp48_lr_b1[j] * x1L + self->hp48_lr_b2[j] * x2L - self->hp48_lr_a1[j] * y1L - self->hp48_lr_a2[j] * y2L;
                    float outR = self->hp48_lr_b0[j] * right_signal + self->hp48_lr_b1[j] * x1R + self->hp48_lr_b2[j] * x2R - self->hp48_lr_a1[j] * y1R - self->hp48_lr_a2[j] * y2R;
                    
                    self->hp48_lr_x2[0][j] = x1L; self->hp48_lr_x1[0][j] = left_signal; self->hp48_lr_y2[0][j] = y1L; self->hp48_lr_y1[0][j] = outL;
                    self->hp48_lr_x2[1][j] = x1R; self->hp48_lr_x1[1][j] = right_signal; self->hp48_lr_y2[1][j] = y1R; self->hp48_lr_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (hp_slope_selector == 5) {
                for (int j = 0; j < 4; ++j) {
                    float x1L = self->hp48_bw_x1[0][j], x2L = self->hp48_bw_x2[0][j], y1L = self->hp48_bw_y1[0][j], y2L = self->hp48_bw_y2[0][j];
                    float x1R = self->hp48_bw_x1[1][j], x2R = self->hp48_bw_x2[1][j], y1R = self->hp48_bw_y1[1][j], y2R = self->hp48_bw_y2[1][j];
                    
                    float outL = self->hp48_bw_b0[j] * left_signal + self->hp48_bw_b1[j] * x1L + self->hp48_bw_b2[j] * x2L - self->hp48_bw_a1[j] * y1L - self->hp48_bw_a2[j] * y2L;
                    float outR = self->hp48_bw_b0[j] * right_signal + self->hp48_bw_b1[j] * x1R + self->hp48_bw_b2[j] * x2R - self->hp48_bw_a1[j] * y1R - self->hp48_bw_a2[j] * y2R;
                    
                    self->hp48_bw_x2[0][j] = x1L; self->hp48_bw_x1[0][j] = left_signal; self->hp48_bw_y2[0][j] = y1L; self->hp48_bw_y1[0][j] = outL;
                    self->hp48_bw_x2[1][j] = x1R; self->hp48_bw_x1[1][j] = right_signal; self->hp48_bw_y2[1][j] = y1R; self->hp48_bw_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (hp_slope_selector == 4) {
                for (int j = 0; j < 2; ++j) {
                    float x1L = self->hp24_bw_x1[0][j], x2L = self->hp24_bw_x2[0][j], y1L = self->hp24_bw_y1[0][j], y2L = self->hp24_bw_y2[0][j];
                    float x1R = self->hp24_bw_x1[1][j], x2R = self->hp24_bw_x2[1][j], y1R = self->hp24_bw_y1[1][j], y2R = self->hp24_bw_y2[1][j];
                    
                    float outL = self->hp24_bw_b0[j] * left_signal + self->hp24_bw_b1[j] * x1L + self->hp24_bw_b2[j] * x2L - self->hp24_bw_a1[j] * y1L - self->hp24_bw_a2[j] * y2L;
                    float outR = self->hp24_bw_b0[j] * right_signal + self->hp24_bw_b1[j] * x1R + self->hp24_bw_b2[j] * x2R - self->hp24_bw_a1[j] * y1R - self->hp24_bw_a2[j] * y2R;
                    
                    self->hp24_bw_x2[0][j] = x1L; self->hp24_bw_x1[0][j] = left_signal; self->hp24_bw_y2[0][j] = y1L; self->hp24_bw_y1[0][j] = outL;
                    self->hp24_bw_x2[1][j] = x1R; self->hp24_bw_x1[1][j] = right_signal; self->hp24_bw_y2[1][j] = y1R; self->hp24_bw_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (hp_slope_selector == 3) {
                for (int j = 0; j < 4; ++j) {
                    float x1L = self->hp48_x1[0][j], x2L = self->hp48_x2[0][j], y1L = self->hp48_y1[0][j], y2L = self->hp48_y2[0][j];
                    float x1R = self->hp48_x1[1][j], x2R = self->hp48_x2[1][j], y1R = self->hp48_y1[1][j], y2R = self->hp48_y2[1][j];
                    
                    float outL = self->hp48_b0[j] * left_signal + self->hp48_b1[j] * x1L + self->hp48_b2[j] * x2L - self->hp48_a1[j] * y1L - self->hp48_a2[j] * y2L;
                    float outR = self->hp48_b0[j] * right_signal + self->hp48_b1[j] * x1R + self->hp48_b2[j] * x2R - self->hp48_a1[j] * y1R - self->hp48_a2[j] * y2R;
                    
                    self->hp48_x2[0][j] = x1L; self->hp48_x1[0][j] = left_signal; self->hp48_y2[0][j] = y1L; self->hp48_y1[0][j] = outL;
                    self->hp48_x2[1][j] = x1R; self->hp48_x1[1][j] = right_signal; self->hp48_y2[1][j] = y1R; self->hp48_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (hp_slope_selector == 2) {
                for (int j = 0; j < 2; ++j) {
                    float x1L = self->hp24_x1[0][j], x2L = self->hp24_x2[0][j], y1L = self->hp24_y1[0][j], y2L = self->hp24_y2[0][j];
                    float x1R = self->hp24_x1[1][j], x2R = self->hp24_x2[1][j], y1R = self->hp24_y1[1][j], y2R = self->hp24_y2[1][j];
                    
                    float outL = self->hp24_b0[j] * left_signal + self->hp24_b1[j] * x1L + self->hp24_b2[j] * x2L - self->hp24_a1[j] * y1L - self->hp24_a2[j] * y2L;
                    float outR = self->hp24_b0[j] * right_signal + self->hp24_b1[j] * x1R + self->hp24_b2[j] * x2R - self->hp24_a1[j] * y1R - self->hp24_a2[j] * y2R;
                    
                    self->hp24_x2[0][j] = x1L; self->hp24_x1[0][j] = left_signal; self->hp24_y2[0][j] = y1L; self->hp24_y1[0][j] = outL;
                    self->hp24_x2[1][j] = x1R; self->hp24_x1[1][j] = right_signal; self->hp24_y2[1][j] = y1R; self->hp24_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else {
                float x1L = self->hp12_x1[0], x2L = self->hp12_x2[0], y1L = self->hp12_y1[0], y2L = self->hp12_y2[0];
                float x1R = self->hp12_x1[1], x2R = self->hp12_x2[1], y1R = self->hp12_y1[1], y2R = self->hp12_y2[1];
                
                float outL = self->hp12_b0 * left_signal + self->hp12_b1 * x1L + self->hp12_b2 * x2L - self->hp12_a1 * y1L - self->hp12_a2 * y2L;
                float outR = self->hp12_b0 * right_signal + self->hp12_b1 * x1R + self->hp12_b2 * x2R - self->hp12_a1 * y1R - self->hp12_a2 * y2R;
                
                self->hp12_x2[0] = x1L; self->hp12_x1[0] = left_signal; self->hp12_y2[0] = y1L; self->hp12_y1[0] = outL;
                self->hp12_x2[1] = x1R; self->hp12_x1[1] = right_signal; self->hp12_y2[1] = y1R; self->hp12_y1[1] = outR;
                
                left_signal = outL; right_signal = outR;
            }
        }
        
        if (lp_slope_selector != 0) {
            if (lp_slope_selector == 8) {
                float x1L = self->hs_x1[0], x2L = self->hs_x2[0], y1L = self->hs_y1[0], y2L = self->hs_y2[0];
                float x1R = self->hs_x1[1], x2R = self->hs_x2[1], y1R = self->hs_y1[1], y2R = self->hs_y2[1];
                
                float outL = self->hs_b0 * left_signal + self->hs_b1 * x1L + self->hs_b2 * x2L - self->hs_a1 * y1L - self->hs_a2 * y2L;
                float outR = self->hs_b0 * right_signal + self->hs_b1 * x1R + self->hs_b2 * x2R - self->hs_a1 * y1R - self->hs_a2 * y2R;
                
                self->hs_x2[0] = x1L; self->hs_x1[0] = left_signal; self->hs_y2[0] = y1L; self->hs_y1[0] = outL;
                self->hs_x2[1] = x1R; self->hs_x1[1] = right_signal; self->hs_y2[1] = y1R; self->hs_y1[1] = outR;
                
                left_signal = outL; right_signal = outR;
            } else if (lp_slope_selector == 7) {
                for (int j = 0; j < 8; ++j) {
                    float x1L = self->lp96_lr_x1[0][j], x2L = self->lp96_lr_x2[0][j], y1L = self->lp96_lr_y1[0][j], y2L = self->lp96_lr_y2[0][j];
                    float x1R = self->lp96_lr_x1[1][j], x2R = self->lp96_lr_x2[1][j], y1R = self->lp96_lr_y1[1][j], y2R = self->lp96_lr_y2[1][j];
                    
                    float outL = self->lp96_lr_b0[j] * left_signal + self->lp96_lr_b1[j] * x1L + self->lp96_lr_b2[j] * x2L - self->lp96_lr_a1[j] * y1L - self->lp96_lr_a2[j] * y2L;
                    float outR = self->lp96_lr_b0[j] * right_signal + self->lp96_lr_b1[j] * x1R + self->lp96_lr_b2[j] * x2R - self->lp96_lr_a1[j] * y1R - self->lp96_lr_a2[j] * y2R;
                    
                    self->lp96_lr_x2[0][j] = x1L; self->lp96_lr_x1[0][j] = left_signal; self->lp96_lr_y2[0][j] = y1L; self->lp96_lr_y1[0][j] = outL;
                    self->lp96_lr_x2[1][j] = x1R; self->lp96_lr_x1[1][j] = right_signal; self->lp96_lr_y2[1][j] = y1R; self->lp96_lr_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (lp_slope_selector == 6) {
                for (int j = 0; j < 4; ++j) {
                    float x1L = self->lp48_lr_x1[0][j], x2L = self->lp48_lr_x2[0][j], y1L = self->lp48_lr_y1[0][j], y2L = self->lp48_lr_y2[0][j];
                    float x1R = self->lp48_lr_x1[1][j], x2R = self->lp48_lr_x2[1][j], y1R = self->lp48_lr_y1[1][j], y2R = self->lp48_lr_y2[1][j];
                    
                    float outL = self->lp48_lr_b0[j] * left_signal + self->lp48_lr_b1[j] * x1L + self->lp48_lr_b2[j] * x2L - self->lp48_lr_a1[j] * y1L - self->lp48_lr_a2[j] * y2L;
                    float outR = self->lp48_lr_b0[j] * right_signal + self->lp48_lr_b1[j] * x1R + self->lp48_lr_b2[j] * x2R - self->lp48_lr_a1[j] * y1R - self->lp48_lr_a2[j] * y2R;
                    
                    self->lp48_lr_x2[0][j] = x1L; self->lp48_lr_x1[0][j] = left_signal; self->lp48_lr_y2[0][j] = y1L; self->lp48_lr_y1[0][j] = outL;
                    self->lp48_lr_x2[1][j] = x1R; self->lp48_lr_x1[1][j] = right_signal; self->lp48_lr_y2[1][j] = y1R; self->lp48_lr_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (lp_slope_selector == 5) {
                for (int j = 0; j < 4; ++j) {
                    float x1L = self->lp48_bw_x1[0][j], x2L = self->lp48_bw_x2[0][j], y1L = self->lp48_bw_y1[0][j], y2L = self->lp48_bw_y2[0][j];
                    float x1R = self->lp48_bw_x1[1][j], x2R = self->lp48_bw_x2[1][j], y1R = self->lp48_bw_y1[1][j], y2R = self->lp48_bw_y2[1][j];
                    
                    float outL = self->lp48_bw_b0[j] * left_signal + self->lp48_bw_b1[j] * x1L + self->lp48_bw_b2[j] * x2L - self->lp48_bw_a1[j] * y1L - self->lp48_bw_a2[j] * y2L;
                    float outR = self->lp48_bw_b0[j] * right_signal + self->lp48_bw_b1[j] * x1R + self->lp48_bw_b2[j] * x2R - self->lp48_bw_a1[j] * y1R - self->lp48_bw_a2[j] * y2R;
                    
                    self->lp48_bw_x2[0][j] = x1L; self->lp48_bw_x1[0][j] = left_signal; self->lp48_bw_y2[0][j] = y1L; self->lp48_bw_y1[0][j] = outL;
                    self->lp48_bw_x2[1][j] = x1R; self->lp48_bw_x1[1][j] = right_signal; self->lp48_bw_y2[1][j] = y1R; self->lp48_bw_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (lp_slope_selector == 4) {
                for (int j = 0; j < 2; ++j) {
                    float x1L = self->lp24_bw_x1[0][j], x2L = self->lp24_bw_x2[0][j], y1L = self->lp24_bw_y1[0][j], y2L = self->lp24_bw_y2[0][j];
                    float x1R = self->lp24_bw_x1[1][j], x2R = self->lp24_bw_x2[1][j], y1R = self->lp24_bw_y1[1][j], y2R = self->lp24_bw_y2[1][j];
                    
                    float outL = self->lp24_bw_b0[j] * left_signal + self->lp24_bw_b1[j] * x1L + self->lp24_bw_b2[j] * x2L - self->lp24_bw_a1[j] * y1L - self->lp24_bw_a2[j] * y2L;
                    float outR = self->lp24_bw_b0[j] * right_signal + self->lp24_bw_b1[j] * x1R + self->lp24_bw_b2[j] * x2R - self->lp24_bw_a1[j] * y1R - self->lp24_bw_a2[j] * y2R;
                    
                    self->lp24_bw_x2[0][j] = x1L; self->lp24_bw_x1[0][j] = left_signal; self->lp24_bw_y2[0][j] = y1L; self->lp24_bw_y1[0][j] = outL;
                    self->lp24_bw_x2[1][j] = x1R; self->lp24_bw_x1[1][j] = right_signal; self->lp24_bw_y2[1][j] = y1R; self->lp24_bw_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (lp_slope_selector == 3) {
                for (int j = 0; j < 4; ++j) {
                    float x1L = self->lp48_x1[0][j], x2L = self->lp48_x2[0][j], y1L = self->lp48_y1[0][j], y2L = self->lp48_y2[0][j];
                    float x1R = self->lp48_x1[1][j], x2R = self->lp48_x2[1][j], y1R = self->lp48_y1[1][j], y2R = self->lp48_y2[1][j];
                    
                    float outL = self->lp48_b0[j] * left_signal + self->lp48_b1[j] * x1L + self->lp48_b2[j] * x2L - self->lp48_a1[j] * y1L - self->lp48_a2[j] * y2L;
                    float outR = self->lp48_b0[j] * right_signal + self->lp48_b1[j] * x1R + self->lp48_b2[j] * x2R - self->lp48_a1[j] * y1R - self->lp48_a2[j] * y2R;
                    
                    self->lp48_x2[0][j] = x1L; self->lp48_x1[0][j] = left_signal; self->lp48_y2[0][j] = y1L; self->lp48_y1[0][j] = outL;
                    self->lp48_x2[1][j] = x1R; self->lp48_x1[1][j] = right_signal; self->lp48_y2[1][j] = y1R; self->lp48_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else if (lp_slope_selector == 2) {
                for (int j = 0; j < 2; ++j) {
                    float x1L = self->lp24_x1[0][j], x2L = self->lp24_x2[0][j], y1L = self->lp24_y1[0][j], y2L = self->lp24_y2[0][j];
                    float x1R = self->lp24_x1[1][j], x2R = self->lp24_x2[1][j], y1R = self->lp24_y1[1][j], y2R = self->lp24_y2[1][j];
                    
                    float outL = self->lp24_b0[j] * left_signal + self->lp24_b1[j] * x1L + self->lp24_b2[j] * x2L - self->lp24_a1[j] * y1L - self->lp24_a2[j] * y2L;
                    float outR = self->lp24_b0[j] * right_signal + self->lp24_b1[j] * x1R + self->lp24_b2[j] * x2R - self->lp24_a1[j] * y1R - self->lp24_a2[j] * y2R;
                    
                    self->lp24_x2[0][j] = x1L; self->lp24_x1[0][j] = left_signal; self->lp24_y2[0][j] = y1L; self->lp24_y1[0][j] = outL;
                    self->lp24_x2[1][j] = x1R; self->lp24_x1[1][j] = right_signal; self->lp24_y2[1][j] = y1R; self->lp24_y1[1][j] = outR;
                    
                    left_signal = outL; right_signal = outR;
                }
            } else {
                float x1L = self->lp12_x1[0], x2L = self->lp12_x2[0], y1L = self->lp12_y1[0], y2L = self->lp12_y2[0];
                float x1R = self->lp12_x1[1], x2R = self->lp12_x2[1], y1R = self->lp12_y1[1], y2R = self->lp12_y2[1];
                
                float outL = self->lp12_b0 * left_signal + self->lp12_b1 * x1L + self->lp12_b2 * x2L - self->lp12_a1 * y1L - self->lp12_a2 * y2L;
                float outR = self->lp12_b0 * right_signal + self->lp12_b1 * x1R + self->lp12_b2 * x2R - self->lp12_a1 * y1R - self->lp12_a2 * y2R;
                
                self->lp12_x2[0] = x1L; self->lp12_x1[0] = left_signal; self->lp12_y2[0] = y1L; self->lp12_y1[0] = outL;
                self->lp12_x2[1] = x1R; self->lp12_x1[1] = right_signal; self->lp12_y2[1] = y1R; self->lp12_y1[1] = outR;
                
                left_signal = outL; right_signal = outR;
            }
        }
        
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
