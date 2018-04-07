#include "index.h"

// class BaseIndex

BaseIndex::BaseIndex(const unsigned index, const IndexMapper* indexMapper) :
  index_(index), indexMapper_(indexMapper) {}

const IndexMapper* BaseIndex::getIndexMapper() const {
  return indexMapper_;
}

unsigned BaseIndex::getIndex() const {
  return index_;
}

// end class BaseIndex

// class BaseIndexWithBlockType

BaseIndexWithBlockType::BaseIndexWithBlockType(const unsigned index, const unsigned blockTypeIndex,
                                               const IndexMapper* indexMapper) :
  BaseIndex(index, indexMapper), blockTypeIndex_(blockTypeIndex) {}

const BlockTypeIndex* BaseIndexWithBlockType::getBlockTypeIndex() const {
  return getIndexMapper()->getBlockTypeIndex(blockTypeIndex_);
}

bool BaseIndexWithBlockType::isDiagonalBlockType() const {
  return getBlockTypeIndex()->isDiagonal();
}

bool BaseIndexWithBlockType::isNonDiagonalBlockType() const {
  return getBlockTypeIndex()->isNonDiagonal();
}

unsigned BaseIndexWithBlockType::getBlockTypeDimension() const {
  return getBlockTypeIndex()->dim_;
}

bool BaseIndexWithBlockType::isBlockTypeDropped() const {
  return getBlockTypeIndex()->isDropped();
}

// end class BaseIndexWithBlockType

// class RandomEffectIndex

RandomEffectIndex::RandomEffectIndex(const unsigned randomEffectIndex,
                                     const unsigned blockTypeIndex,
                                     const IndexMapper* indexMapper,
                                     const unsigned blockIndex) :
  BaseIndexWithBlockType(randomEffectIndex, blockTypeIndex, indexMapper),
  blockIndex_(blockIndex) {}

const BlockIndex* RandomEffectIndex::getBlockIndex() const {
  return getIndexMapper()->getBlockIndex(blockIndex_);
}

// end class RandomEffectIndex

// class BlockIndex

BlockIndex::BlockIndex(const unsigned blockIndex, const unsigned blockTypeIndex,
                       const IndexMapper* indexMapper,
                       const std::vector<RandomEffectIndex*> randomEffects) :
  BaseIndexWithBlockType(blockIndex, blockTypeIndex, indexMapper),
  randomEffects_(randomEffects) {}

const std::vector<RandomEffectIndex*>& BlockIndex::getRandomEffects() const {
  return randomEffects_;
}

// end class BlockIndex

// class ThetaIndex

ThetaIndex::ThetaIndex(const unsigned thetaIndex, const unsigned blockTypeIndex,
                       const IndexMapper* indexMapper) :
  BaseIndexWithBlockType(thetaIndex, blockTypeIndex, indexMapper) {}

// end class ThetaIndex

// class BlockTypeIndex

BlockTypeIndex::BlockTypeIndex(const unsigned blockTypeIndex, const IndexMapper* indexMapper,
                               const unsigned dim,
                               const std::vector<RandomEffectIndex*> randomEffects,
                               const std::vector<BlockIndex*> blocks,
                               const std::vector<ThetaIndex*> thetas) :
  BaseIndex(blockTypeIndex, indexMapper), dim_(dim), blockTypeDropped_(false),
  randomEffects_(randomEffects), blocks_(blocks), thetas_(thetas) {}

bool BlockTypeIndex::isNonDiagonal() const {
  return dim_ > 1;
}

bool BlockTypeIndex::isDiagonal() const {
  return dim_ == 1;
}

bool BlockTypeIndex::isDropped() const {
  return blockTypeDropped_;
}

const std::vector<RandomEffectIndex*>& BlockTypeIndex::getRandomEffects() const {
  return randomEffects_;
}

const std::vector<BlockIndex*>& BlockTypeIndex::getBlocks() const {
  return blocks_;
}

const std::vector<ThetaIndex*>& BlockTypeIndex::getThetas() const {
  return thetas_;
}

unsigned BlockTypeIndex::getNumberOfBlocks() const {
  return blocks_.size();
}

bool BlockTypeIndex::isDropped(const VectorXd& theta) const {
  for (const ThetaIndex* thetaIndex : thetas_) {
    if (theta[thetaIndex->getIndex()] != 0.) { return false; }
  }
  return true;
}

void BlockTypeIndex::updateDropped(const VectorXd& theta) {
  blockTypeDropped_ = isDropped(theta);
}

// end class BlockTypeIndex

// class IndexMapper

IndexMapper::IndexMapper(const Rcpp::List& args) :
  ind_(Rcpp::as<MiVec>(args["ind"])), k_(Rcpp::as<MiVec>(args["k"])), bBlockMap_(repDo(ind_, k_)),
  anyBlockTypeDropped_(false), allBlockTypesDropped_(false) {
  MiVec dim(Rcpp::as<MiVec>(args["dim"]));
  std::vector<VectorXi> idx(Rcpp::as<std::vector<VectorXi> >(Rcpp::as<Rcpp::List>(args["idx"]))),
  blockBMap(Rcpp::as<std::vector<VectorXi> >(Rcpp::as<Rcpp::List>(args["blockBMap"]))),
  blockThetaMap(Rcpp::as<std::vector<VectorXi> >(Rcpp::as<Rcpp::List>(args["thetaBlockMap"])));

  initialiseRandomEffects();
  initialiseBlocks(blockBMap);
  initialiseThetas(blockThetaMap);
  initialiseBlockTypes(dim, idx, blockThetaMap);
}

void IndexMapper::initialiseRandomEffects() {
  randomEffects_.reserve(k_.size());
  for (int randomEffect = 0, size = k_.size(); randomEffect < size; ++randomEffect) {
    int block(k_[randomEffect] - 1);
    int blockType(ind_[block] - 1);
    randomEffects_.emplace_back(new RandomEffectIndex(randomEffect, blockType, this, block));
  }
}

void IndexMapper::initialiseBlocks(const std::vector<VectorXi>& blockBMap) {
  blocks_.reserve(ind_.size());
  for (int block = 0, size = ind_.size(); block < size; ++block) {
    int blockTypeIndex(ind_[block] - 1);
    std::vector<RandomEffectIndex*> randomEffects;
    randomEffects.reserve(blockBMap[block].size());
    for (int randomEffect = 0, blockSize = blockBMap[block].size();
         randomEffect < blockSize; ++randomEffect) {
      randomEffects.push_back(randomEffects_[blockBMap[block][randomEffect] - 1].get());
    }
    blocks_.emplace_back(new BlockIndex(block, blockTypeIndex, this, randomEffects));
  }
}

void IndexMapper::initialiseThetas(const std::vector<VectorXi>& blockThetaMap) {
  int nThetas = 0;
  for (int blockType = 0, size = blockThetaMap.size(); blockType < size; ++blockType) {
    nThetas +=  blockThetaMap[blockType].size();
  }
  thetas_.reserve(nThetas);
  for (int blockType = 0, size = blockThetaMap.size(); blockType < size; ++blockType) {
    for (int theta = 0, nThetasForBlock = blockThetaMap[blockType].size(); theta < nThetasForBlock; ++theta) {
      thetas_.emplace_back(new ThetaIndex(blockThetaMap[blockType][theta] - 1, blockType, this));
    }
  }
}

void IndexMapper::initialiseBlockTypes(const MiVec& dim, const std::vector<VectorXi>& idx,
                                       const std::vector<VectorXi>& blockThetaMap) {
  anyBlockTypeNonDiagonal_ = false;
  blockTypes_.reserve(dim.size());
  for (int blockType = 0, size = dim.size(); blockType < size; ++blockType) {
    if (dim[blockType] > 1) {
      anyBlockTypeNonDiagonal_ = true;
    }

    std::vector<RandomEffectIndex*> randomEffects;
    randomEffects.reserve(idx[blockType].size());
    for (int randomEffect = 0, nRandomEffects = idx[blockType].size();
         randomEffect < nRandomEffects; ++randomEffect) {
      randomEffects.push_back(randomEffects_[idx[blockType][randomEffect] - 1].get());
    }

    std::vector<BlockIndex*> blocks;
    blocks.reserve(randomEffects.size() / dim[blockType]);
    for (int block = 0, nBlocks = ind_.size(); block < nBlocks; ++block) {
      if (ind_[block] == blockType + 1) {
        blocks.push_back(blocks_[block].get());
      }
    }

    std::vector<ThetaIndex*> thetas;
    thetas.reserve(blockThetaMap[blockType].size());
    for (int theta = 0, nThetas = blockThetaMap[blockType].size(); theta < nThetas; ++theta) {
      thetas.push_back(thetas_[blockThetaMap[blockType][theta] - 1].get());
    }

    blockTypes_.emplace_back(new BlockTypeIndex(blockType, this, dim[blockType], randomEffects,
                                                blocks, thetas));
  }

}

const MiVec& IndexMapper::getInd() const {
  return ind_;
}

const MiVec& IndexMapper::getK() const {
  return k_;
}

const VectorXi& IndexMapper::getBBlockMap() const {
  return bBlockMap_;
}

unsigned IndexMapper::getNumberOfRandomEffects() const {
  return randomEffects_.size();
}

unsigned IndexMapper::getNumberOfBlocks() const {
  return blocks_.size();
}

unsigned IndexMapper::getNumberOfThetas() const {
  return thetas_.size();
}

unsigned IndexMapper::getNumberOfBlockTypes() const {
  return blockTypes_.size();
}

const std::vector<std::unique_ptr<RandomEffectIndex> >& IndexMapper::getRandomEffects() const {
  return randomEffects_;
}

const std::vector<std::unique_ptr<BlockIndex> >& IndexMapper::getBlocks() const {
  return blocks_;
}

const std::vector<std::unique_ptr<ThetaIndex> >& IndexMapper::getThetas() const {
  return thetas_;
}

const std::vector<std::unique_ptr<BlockTypeIndex> >& IndexMapper::getBlockTypes() const {
  return blockTypes_;
}

std::vector<std::unique_ptr<BlockTypeIndex> >& IndexMapper::getBlockTypes() {
  return blockTypes_;
}

unsigned IndexMapper::getMaxBlockTypeDimension() const {
  unsigned dim = 0;
  for (const auto &blockType : blockTypes_) {
    if (blockType->dim_ > dim) dim = blockType->dim_;
  }
  return dim;
}

bool IndexMapper::isAnyBlockTypeNonDiagonal() const {
  return anyBlockTypeNonDiagonal_;
}

bool IndexMapper::isAnyBlockTypeDropped() const {
  return anyBlockTypeDropped_;
}

bool IndexMapper::areAllBlockTypesDropped() const {
  return allBlockTypesDropped_;
}

const RandomEffectIndex* IndexMapper::getRandomEffect(unsigned randomEffectIndex) const {
  return randomEffects_[randomEffectIndex].get();
}

const BlockIndex* IndexMapper::getBlockIndex(unsigned blockIndex) const {
  return blocks_[blockIndex].get();
}

const ThetaIndex* IndexMapper::getThetaIndex(unsigned thetaIndex) const {
  return thetas_[thetaIndex].get();
}

const BlockTypeIndex* IndexMapper::getBlockTypeIndex(unsigned blockTypeIndex) const {
  return blockTypes_[blockTypeIndex].get();
}

void IndexMapper::updateBlockTypesDropped(const VectorXd& theta) {
  allBlockTypesDropped_ = true, anyBlockTypeDropped_ = false;
  for (auto &blockType : blockTypes_) {
    blockType->updateDropped(theta);
    if (!blockType->isDropped()) {
      allBlockTypesDropped_ = false;
    } else {
      anyBlockTypeDropped_ = true;
    }
  }
}

// end class IndexMapper
