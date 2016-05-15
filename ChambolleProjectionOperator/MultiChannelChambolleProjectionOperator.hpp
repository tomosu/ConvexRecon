#pragma once

std::vector<Reconstruction> MultiChannelChambolleProjectionOperator(const Real_t lambda, const Real_t step, const Real_t tolerance, const int maxIteration,
								    const Geometry &geo, const std::vector<Reconstruction> &input);
