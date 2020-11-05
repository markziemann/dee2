module ResultsPage.Types exposing (..)

import Array
import Dict
import SearchPage.Types
import SharedTypes
import Table
import Set

type alias SelectedResults =
    Dict.Dict Int ( String, String )


type alias ResultsPendingRemoval =
    Set.Set Int


type alias Model =
    { searchHits : Maybe Int
    , searchResultRows : Maybe (Array.Array SearchPage.Types.SearchResult)
    , resultsTableQuery : String
    , resultsTableState : Table.State
    , selectedResultsTableQuery : String
    , selectedResultsTableState : Table.State
    , downloading : Bool
    , selectedResults : SelectedResults
    , resultsPendingRemoval : ResultsPendingRemoval
    , paginationOffset : SharedTypes.PaginationOffset
    }


type Msg
    = ResultClicked SearchPage.Types.SearchResult
    | RemoveStagedSelections
    | SelectedResultClicked Int
    | SetResultsTableQuery String
    | SetResultsTableState Table.State
    | SetSelectedResultsTableQuery String
    | SetSelectedResultsTableState Table.State
    | DownloadRequested
    | DownloadButtonReset


type OutMsg
    = RequestNextPage SharedTypes.PaginationOffset
