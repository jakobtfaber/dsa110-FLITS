import numpy as np

from scattering.scat_analysis.pipeline.io import BurstDataset


def dataset_for_crop(data):
    dataset = BurstDataset.__new__(BurstDataset)
    dataset.name = "test"
    dataset.data = np.asarray(data, dtype=float)
    dataset.time = np.arange(dataset.data.shape[1], dtype=float)
    dataset.dt_ms = 1.0
    dataset.smooth_ms = 0.0
    dataset.onpulse_thresh = 3.0
    dataset.onpulse_pad_factor = 0.5
    dataset.onpulse_crop_status = "not_requested"
    dataset.onpulse_crop_bounds = None
    return dataset


def test_crop_status_fails_closed_when_no_sample_clears_threshold():
    dataset = dataset_for_crop(np.zeros((4, 16)))
    dataset._crop_on_pulse()
    assert dataset.onpulse_crop_status == "failed_zero_offpulse_spread"
    assert dataset.onpulse_crop_bounds is None


def test_crop_status_records_applied_bounds():
    data = np.zeros((4, 32))
    data[:, 15:17] = 20.0
    # Give the off-pulse region non-zero robust spread.
    data[:, :8] = np.tile([-1.0, 1.0], 4)
    data[:, -8:] = np.tile([1.0, -1.0], 4)
    dataset = dataset_for_crop(data)
    dataset._crop_on_pulse()
    assert dataset.onpulse_crop_status == "applied"
    assert dataset.onpulse_crop_bounds is not None
    assert dataset.data.shape[1] < data.shape[1]


def test_crop_status_fails_closed_for_single_sample_crop():
    data = np.zeros((4, 32))
    data[:, 16] = 20.0
    data[:, :8] = np.tile([-1.0, 1.0], 4)
    data[:, -8:] = np.tile([1.0, -1.0], 4)
    dataset = dataset_for_crop(data)
    dataset.onpulse_pad_factor = 0.25
    dataset._crop_on_pulse()
    assert dataset.onpulse_crop_status == "failed_insufficient_crop_samples"
    assert dataset.onpulse_crop_bounds is None
    assert dataset.data.shape == data.shape
